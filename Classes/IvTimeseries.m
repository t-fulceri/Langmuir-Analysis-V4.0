classdef IvTimeseries < handle
    %An object of this class represent a timeseries of IV langmuir curves
    %and it is created starting from an IvSignals object
    
    properties
        ivs_source_plasma
        ivs_source_vacuum
        V_cycles
        I_cycles
        NC_tot
        time_steps_per_cycle
        initial_time_steps_offset
        V_cycles_rise
        V_cycles_fall
        I_cycles_rise
        I_cycles_fall
        time_uncertainty
        total_time_steps
        total_time_axis_rise
        total_time_axis_fall
        rise
        fall
    end
    
    methods
        function obj = IvTimeseries(varargin)
            %Expecting one of the following:
            %IvTimeseries(plasma_full_file_path)
            %IvTimeseries(plasma_full_file_path,vacuum_full_file_path)
            switch nargin
                case 0
                    %No input arguments -> error
                    fprintf('Error: No input data!\n');
                    return
                case 1
                    obj.ivs_source_plasma = IvSignals(varargin{1});
                    obj.ivs_source_vacuum = IvSignals;
                case 2
                    obj.ivs_source_plasma = IvSignals(varargin{1});
                    obj.ivs_source_vacuum = IvSignals(varargin{2});
            end
            
            %obj.ivs_source_plasma.calculateCMEffect;
            %Split the V,I signals into chunks of equal duration
            %(each corresponding to one voltage cycle)
            voltage = obj.ivs_source_plasma.V_signal;
            
%             %Optional removal of vacuum effects
%             if obj.ivs_source_vacuum.init
%                 current = obj.ivs_source_plasma.I_CM_corrected - obj.ivs_source_vacuum.extractVacuumEffect;
%             else
%                 current = obj.ivs_source_plasma.I_CM_corrected;
%             end
            %Optional removal of vacuum effects
            if obj.ivs_source_vacuum.init
                current = obj.ivs_source_plasma.I_signal - obj.ivs_source_vacuum.extractVacuumEffect;
            else
                current = obj.ivs_source_plasma.I_signal;
            end

            
            f_sampling = obj.ivs_source_plasma.f_sampling;
            time_step = obj.ivs_source_plasma.time_step;
            [obj.V_cycles,obj.initial_time_steps_offset] =  timeseriesSplitSinefit(voltage,voltage,f_sampling);
            obj.I_cycles =  timeseriesSplitSinefit(voltage,current,f_sampling);

            %Total number of cycles
            obj.NC_tot = size(obj.V_cycles,2);
            %Time steps per cycle
            obj.time_steps_per_cycle = size(obj.V_cycles,1);

            hl = round(obj.time_steps_per_cycle/2);

            obj.V_cycles_rise = obj.V_cycles(1:hl,:);
            obj.I_cycles_rise = obj.I_cycles(1:hl,:);

            obj.V_cycles_fall = obj.V_cycles(hl+1:end,:);
            obj.I_cycles_fall = obj.I_cycles(hl+1:end,:);


            %obj.time_steps_per_cycle = time_steps_per_cycle_array(1);
            obj.time_uncertainty = round(obj.time_steps_per_cycle/4)*time_step;

            %obj.initial_time_steps_offset = obj.initial_time_steps_offset_array(1);


            %Total number of time steps
            obj.total_time_steps = obj.initial_time_steps_offset + obj.time_steps_per_cycle*obj.NC_tot;
            fprintf('Total number of time steps extracted\n');
            %Total time axis for rising edge
            obj.total_time_axis_rise = obj.ivs_source_plasma.time_step*(obj.initial_time_steps_offset + round(obj.time_steps_per_cycle/4) + (0:obj.NC_tot-1)*obj.time_steps_per_cycle);
            fprintf('Total time axis for rising edge extracted\n');
            %Total time axis for falling edge
            obj.total_time_axis_fall = obj.ivs_source_plasma.time_step*(obj.initial_time_steps_offset + round(obj.time_steps_per_cycle*(3/4)) + (0:obj.NC_tot-1)*obj.time_steps_per_cycle);
            fprintf('Total time axis for falling edge extracted\n');

            %Stack together the datapoints and data uncertainties which correspond to
            %the same cycle
            iv_cycles_stack = struct('V_cycles_rise',[],'I_cycles_rise',[],'V_cycles_fall',[],'I_cycles_fall',[]);

            iv_cycles_stack.V_cycles_rise = cat(1,iv_cycles_stack.V_cycles_rise,obj.V_cycles_rise);
            iv_cycles_stack.I_cycles_rise = cat(1,iv_cycles_stack.I_cycles_rise,obj.I_cycles_rise);
            iv_cycles_stack.V_cycles_fall = cat(1,iv_cycles_stack.V_cycles_fall,obj.V_cycles_fall);
            iv_cycles_stack.I_cycles_fall = cat(1,iv_cycles_stack.I_cycles_fall,obj.I_cycles_fall);

            fprintf('Datapoints have been stacked together\n');

            %Sort the V(t) and I(t) signals in every cycle following the V-order
            for nc = 1:obj.NC_tot
                [iv_cycles_stack.V_cycles_rise(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_rise(:,nc));
                iv_cycles_stack.I_cycles_rise(:,nc) = iv_cycles_stack.I_cycles_rise(sort_ind,nc);

                [iv_cycles_stack.V_cycles_fall(:,nc),sort_ind] = sort(iv_cycles_stack.V_cycles_fall(:,nc));
                iv_cycles_stack.I_cycles_fall(:,nc) = iv_cycles_stack.I_cycles_fall(sort_ind,nc);
            end
            fprintf('V(t) and I(t) signals sorted according to V-order\n');


            for nc = 1:obj.NC_tot

                rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_rise(:,nc),iv_cycles_stack.I_cycles_rise(:,nc));
                obj.rise(nc).V_axis = rec_curve_out.x;
                obj.rise(nc).I_of_V = rec_curve_out.y_of_x;

                rec_curve_out = reconstructCurve(iv_cycles_stack.V_cycles_fall(:,nc),iv_cycles_stack.I_cycles_fall(:,nc));
                obj.fall(nc).V_axis = rec_curve_out.x;
                obj.fall(nc).I_of_V = rec_curve_out.y_of_x;

            end
            fprintf('IV-timeseries ready\n');

            end
    end
    
end

