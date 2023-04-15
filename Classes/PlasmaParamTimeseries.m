classdef PlasmaParamTimeseries < handle
    %PLASMAPARAMTIMESERIES Contains the information on the time evolution of plasma
    %parameters reconstructed by fitting a timeseries of IV curves with a
    %simple model for le tha Langmuir IV charactersitic
    %Tiziano Fulceri 2019-02-25
    
    properties
        ivt_source
        vft_source
        
        plasma_param_extracted = false
        
        fitresult_ca
        gof_ca
        I_of_V_fitted_ca
        I_ion_fitted_ca
        diagnostics_ca
        
        T_e_a
        T_e_upper_a
        T_e_lower_a
        
        V_float_a
        V_plasma_a
        
        I_sat_e_a
        I_sat_e_upper_a
        I_sat_e_lower_a
        
        I_sat_i_a
        I_sat_i_upper_a
        I_sat_i_lower_a
    end
    
    methods
        
        function obj = PlasmaParamTimeseries(varargin)
            %Expecting one of the following:
            %PlasmaParamTimeseries(plasma_full_file_path)
            %PlasmaParamTimeseries(plasma_full_file_path,v_float_full_file_path)
            
            if isempty(varargin)
                    %No input arguments -> error
                    fprintf('Error: No input data!\n');
                    return
            end
            
            if nargin >= 1 && exist(varargin{1}) %#ok<EXIST>
                obj.ivt_source = IvTimeseries(varargin{1});
            end
            if nargin >= 2 && exist(varargin{2}) %#ok<EXIST>
                obj.vft_source = VTimeseries(obj.ivt_source,varargin{2});
            end
            
            Nc = obj.ivt_source.NC_tot;
            
            obj.fitresult_ca = cell(Nc,2);
            obj.gof_ca = cell(Nc,2);
            obj.I_of_V_fitted_ca = cell(Nc,2);
            obj.I_ion_fitted_ca = cell(Nc,2);
            obj.diagnostics_ca = cell(Nc,2);
            
            obj.T_e_a = zeros(Nc,2);
            obj.T_e_upper_a = zeros(Nc,2);
            obj.T_e_lower_a = zeros(Nc,2);
            
            obj.V_float_a = zeros(Nc,2);
            obj.V_plasma_a = zeros(Nc,2);
            obj.I_sat_e_a = zeros(Nc,2);
            obj.I_sat_i_a = zeros(Nc,2);
            
            fprintf('PlasmaParameter-timeseries object created\n');
        end
        
        function obj = extractPlasmaParameters(obj,S_probe,A_ion,langmuir_iv_fit_version,options,varargin)
            
            if isfield(options,'diag_msg_flag')
                diag_msg_flag = options.diag_msg_flag;
            else
                diag_msg_flag = false;
            end
            
            switch length(varargin)
                case 0
                    nc_range = 1:obj.ivt_source.NC_tot;
                    rf_range = 1:2;
                case 1
                    nc_range = varargin{1};
                    rf_range = 1:2;
                case 2
                    nc_range = varargin{1};
                    rf_range = varargin{2};
            end
            
            %Preparation loop
            %For every voltage cycle
            for nc = nc_range
                %For every rising and faling edge
                for rf = rf_range
                %Output parameters from rising edges
                    obj.fitresult_ca{nc,rf} = [];
                    obj.gof_ca{nc,rf} = [];
                    obj.I_of_V_fitted_ca{nc,rf} = [];
                    obj.I_ion_fitted_ca{nc,rf} = [];

                    obj.T_e_a(nc,rf) = nan;
                    obj.T_e_upper_a(nc,rf) = nan;
                    obj.T_e_lower_a(nc,rf) = nan;

                    obj.V_float_a(nc,rf) = nan;
                    obj.V_plasma_a(nc,rf) = nan;
                    obj.I_sat_i_a(nc,rf) = nan;
                end
            end
            
            %Main loop
            %For every voltage cycle
            for  nc = nc_range
                if diag_msg_flag
                    fprintf(['Plasma parameter extraction of cycle number ',num2str(nc),' out of ',num2str(obj.ivt_source.NC_tot),' starting...\n']);
                end
                %For every rising and faling edge
                for rf = rf_range

                    switch rf
                        case 1 %Rise
                            if diag_msg_flag
                                fprintf('(Rising edge)\n');
                            end
                            V_axis = obj.ivt_source.rise(nc).V_axis;
                            I_of_V = obj.ivt_source.rise(nc).I_of_V;
                        case 2 %Fall
                            if diag_msg_flag
                                fprintf('(Falling edge)\n');
                            end
                            V_axis = obj.ivt_source.fall(nc).V_axis;
                            I_of_V = obj.ivt_source.fall(nc).I_of_V;
                    end
                    
                    %Execute the Langmuir IV fit on the selected IV arrays
                    [return_status, fit_iv_curve_out, diagnostics] = LangmuirIvFit.fitIvCurve(V_axis,I_of_V,S_probe,A_ion,langmuir_iv_fit_version,options);
                    obj.diagnostics_ca{nc,rf} = diagnostics;
                    
                    if return_status == 1
                        
                        obj.fitresult_ca{nc,rf} = fit_iv_curve_out.fitresult;
                        obj.gof_ca{nc,rf} = fit_iv_curve_out.gof;
                        obj.I_of_V_fitted_ca{nc,rf} = fit_iv_curve_out.I_of_V_fitted;
                        obj.I_ion_fitted_ca{nc,rf} = fit_iv_curve_out.I_ion_fitted;
                        

                        obj.T_e_a(nc,rf) = fit_iv_curve_out.T_e;
                        obj.T_e_upper_a(nc,rf) = fit_iv_curve_out.T_e_upper;
                        obj.T_e_lower_a(nc,rf) = fit_iv_curve_out.T_e_lower;

                        obj.V_float_a(nc,rf) = fit_iv_curve_out.V_float;
                        obj.V_plasma_a(nc,rf) = fit_iv_curve_out.V_plasma;
                        
                        obj.I_sat_e_a(nc,rf) = fit_iv_curve_out.I_sat_e;
                        obj.I_sat_e_upper_a(nc,rf) = fit_iv_curve_out.I_sat_e_upper;
                        obj.I_sat_e_lower_a(nc,rf) = fit_iv_curve_out.I_sat_e_lower;

                        obj.I_sat_i_a(nc,rf) = fit_iv_curve_out.I_sat_i;
                        obj.I_sat_i_upper_a(nc,rf) = fit_iv_curve_out.I_sat_i_upper;
                        obj.I_sat_i_lower_a(nc,rf) = fit_iv_curve_out.I_sat_i_lower;
                    end
                end
                if diag_msg_flag
                    fprintf(['End of plasma parameter extraction of cycle number ',num2str(nc),' out of ',num2str(obj.ivt_source.NC_tot),'\n\n']);
                end
            end 
            
        end
        
    end
    
end

