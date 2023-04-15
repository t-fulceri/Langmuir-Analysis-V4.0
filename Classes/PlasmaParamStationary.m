classdef PlasmaParamStationary < handle
    %PLASMAPARAMSTATIONARY Contains the information on sationary plasma parameters
    %extracted from a time-averaged IvTimeseries. It is designed to be used
    %for the analysis os stationary plasma measurement (such as the ECRH)
    %Tiziano Fulceri 2019-03-06
    
    properties
        ivt_source
        vst_source
        
        plasma_param_extracted = false
        
        V_axis_rise;
        I_of_V_rise;
        I_of_V_rise_STD;
        
        V_float_rise;
        
        V_axis_fall;
        I_of_V_fall;
        I_of_V_fall_STD;
        
        V_float_fall;
        
        fit_iv_curve_out_ca = cell(1,2);
    end
    
    methods
        
        function obj = PlasmaParamStationary(varargin)
            %Expecting one of the following:
            %PlasmaParamStationary(iv_signals_full_file_path,v_float_signal_full_file_path)
            
            if isempty(varargin)
                %No input arguments -> error
                fprintf('Error: No input data!\n');
                return
            end
            
            if exist(varargin{1}) %#ok<EXIST>
                obj.ivt_source = IvTimeseries(varargin{1});
                obj.buildViRiseFall;
            else
                fprintf('No IV timeseries input!');
                return
            end
            
            
            if exist(varargin{2}) %#ok<EXIST>
                obj.vst_source = VTimeseries(varargin{1},varargin{2});
                obj.buildVFloatRiseFall;
            else
                fprintf('No V_float timeseries input: uninitialized parameter');
            end
            
            fprintf('PlasmaParamStationary object created\n');
        end
        
        function obj = extractPlasmaParameters(obj,S_probe,A_ion,langmuir_iv_fit_version,options)
            for rf = 1:2
                switch rf
                    case 1
                        options.V_bars = (obj.ivt_source.ivs_source_plasma.V_ADC_step/2).*ones(length(obj.V_axis_rise),1);
                        options.I_bars = obj.I_of_V_rise_STD;
                        if options.include_external_V_float
                            options.external_V_float = obj.V_float_rise;
                        end
                        [return_status, fit_iv_curve_out, ~] = LangmuirIvFit.fitIvCurve(obj.V_axis_rise,obj.I_of_V_rise,S_probe,A_ion,langmuir_iv_fit_version,options);
                    case 2
                        options.V_bars = (obj.ivt_source.ivs_source_plasma.V_ADC_step/2).*ones(length(obj.V_axis_fall),1);
                        options.I_bars = obj.I_of_V_fall_STD;
                        if options.include_external_V_float
                            options.external_V_float = obj.V_float_fall;
                        end
                        [return_status, fit_iv_curve_out, ~] = LangmuirIvFit.fitIvCurve(obj.V_axis_fall,obj.I_of_V_fall,S_probe,A_ion,langmuir_iv_fit_version,options);
                end
                if return_status == 1
                    obj.fit_iv_curve_out_ca{rf} = fit_iv_curve_out;
                end
            end
            
        end
        
        function obj = buildViRiseFall(obj)
            
            V_cycles_rise = obj.ivt_source.V_cycles_rise;
            I_cycles_rise = obj.ivt_source.I_cycles_rise;
            s = size(V_cycles_rise);
            
            V_values_rise = reshape(V_cycles_rise,1,s(1)*s(2));
            I_values_rise = reshape(I_cycles_rise,1,s(1)*s(2));
            
            V_cycles_fall = obj.ivt_source.V_cycles_fall;
            I_cycles_fall = obj.ivt_source.I_cycles_fall;
            s = size(V_cycles_fall);
            
            V_values_fall = reshape(V_cycles_fall,1,s(1)*s(2));
            I_values_fall = reshape(I_cycles_fall,1,s(1)*s(2));
            
            
            rec_curve_out = reconstructCurve(V_values_rise,I_values_rise);
            obj.V_axis_rise = rec_curve_out.x;
            obj.I_of_V_rise = rec_curve_out.y_of_x;
            obj.I_of_V_rise_STD = rec_curve_out.std_y;
            
            rec_curve_out = reconstructCurve(V_values_fall,I_values_fall);
            obj.V_axis_fall = rec_curve_out.x;
            obj.I_of_V_fall = rec_curve_out.y_of_x;
            obj.I_of_V_fall_STD = rec_curve_out.std_y;
            
            
            
        end
        
        function obj = buildVFloatRiseFall(obj)
            
            V_cycles_rise = obj.vst_source.V_cycles_rise;
            s = size(V_cycles_rise);
            obj.V_float_rise = reshape(V_cycles_rise,1,s(1)*s(2));
            
            V_cycles_fall = obj.vst_source.V_cycles_fall;
            s = size(V_cycles_fall);
            obj.V_float_fall = reshape(V_cycles_fall,1,s(1)*s(2));
            
        end
    end
    
end

