classdef VTimeseries < handle
    %An object of this class represent a timeseries of IV langmuir curves
    %and it is created starting from an IvSignals object
    
    properties
        ivs_time_reference_source
        vs_source
        V_cycles
        NC_tot
        time_steps_per_cycle
        initial_time_steps_offset
        V_cycles_rise
        V_cycles_fall
        time_uncertainty
        total_time_steps
        total_time_axis_rise
        total_time_axis_fall
        rise
        fall
    end
    
    methods
        function obj = VTimeseries(varargin)
            %VTimeseries(vs_source,ivs_time_reference_source)
            if isempty(varargin)
                %No input arguments -> error
                fprintf('Error: No input data!\n');
                return
            end
            
            if exist(varargin{1}) %#ok<EXIST>
                obj.ivs_time_reference_source = IvSignals(varargin{1});
            else
                printf('No IV signals reference source folder!');
            end
            
            if exist(varargin{2}) %#ok<EXIST>
                obj.vs_source = VSignal(varargin{2});
            else
                printf('No voltage signal source folder!');
            end
            

            
            voltage = obj.vs_source.V_signal;
            time_reference = obj.ivs_time_reference_source.V_signal;
 
            f_sampling = obj.ivs_time_reference_source.f_sampling;
            time_step = obj.ivs_time_reference_source.time_step;
            
            [obj.V_cycles,obj.initial_time_steps_offset] =  timeseriesSplitSinefit(time_reference,voltage,f_sampling);

            %Total number of cycles
            obj.NC_tot = size(obj.V_cycles,2);
            %Time steps per cycle
            obj.time_steps_per_cycle = size(obj.V_cycles,1);

            hl = round(obj.time_steps_per_cycle/2);

            obj.V_cycles_rise = obj.V_cycles(1:hl,:);

            obj.V_cycles_fall = obj.V_cycles(hl+1:end,:);

            %obj.time_steps_per_cycle = time_steps_per_cycle_array(1);
            obj.time_uncertainty = round(obj.time_steps_per_cycle/4)*time_step;

            %obj.initial_time_steps_offset = obj.initial_time_steps_offset_array(1);

            %Total number of time steps
            obj.total_time_steps = obj.initial_time_steps_offset + obj.time_steps_per_cycle*obj.NC_tot;
            fprintf('Total number of time steps extracted\n');
            %Total time axis for rising edge
            obj.total_time_axis_rise = obj.ivs_time_reference_source.time_step*(obj.initial_time_steps_offset + round(obj.time_steps_per_cycle/4) + (0:obj.NC_tot-1)*obj.time_steps_per_cycle);
            fprintf('Total time axis for rising edge extracted\n');
            %Total time axis for falling edge
            obj.total_time_axis_fall = obj.ivs_time_reference_source.time_step*(obj.initial_time_steps_offset + round(obj.time_steps_per_cycle*(3/4)) + (0:obj.NC_tot-1)*obj.time_steps_per_cycle);
            fprintf('Total time axis for falling edge extracted\n');

            fprintf('V-timeseries ready\n');

            end
    end
    
end

