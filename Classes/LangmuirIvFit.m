classdef LangmuirIvFit < handle
    %LANGMUIRIVFIT Contains all the possible fitting procedures for a
    %Langmuir IV characteristic
    
    properties (Constant)
    end
    
    methods (Static)
        
        %GENERAL CALL - SELECT VERSION BY GIVING IT AS INPUT ARGUMENT
        function [return_status, fit_iv_curve_out, diagnostics] = fitIvCurve(V_axis,I_of_V,S_probe,A_ion,version,options)
            
            switch version
                case 1
                    [return_status, fit_iv_curve_out, diagnostics] = LangmuirIvFit.fitIvCurveVers01(V_axis,I_of_V,S_probe,A_ion,options);
                case 2
                    [return_status, fit_iv_curve_out, diagnostics] = LangmuirIvFit.fitIvCurveVers02(V_axis,I_of_V,S_probe,A_ion,options);
                case 3
                    [return_status, fit_iv_curve_out, diagnostics] = LangmuirIvFit.fitIvCurveVers03(V_axis,I_of_V,S_probe,A_ion,options);
            end
            
        end

        %% VERSION 1 - STABLE - GOOD FOR HIGH DENSITY PLASMA
        %Issue: V_plasma is offsetted to the right, I_sat_i, T_e, V_plasma,and V_float are not mutually consistent (Olaf's comment 2019-03-11)
        %(n_e ~10^{18}%10^{19} m{^-3})
        function [return_status, fit_iv_curve_out, diagnostics] = fitIvCurveVers01(V_axis,I_of_V,S_probe,A_ion,options)

            fprintf('Running fitIvCurveVers01 with following parameters:\n\n');
            %Prepare return status flag
            return_status = 0;
            %Initialize the fit_iv_curve_out output structure
            fit_iv_curve_out.fitresult_cell_array = [];
            fit_iv_curve_out.gof_cell_array = [];
            fit_iv_curve_out.I_of_V_fitted = [];
            fit_iv_curve_out.I_ion_fitted = [];

            fit_iv_curve_out.T_e = nan;
            fit_iv_curve_out.T_e_upper = nan;
            fit_iv_curve_out.T_e_lower = nan;

            fit_iv_curve_out.V_float = nan;
            fit_iv_curve_out.V_plasma = nan;
            
            fit_iv_curve_out.I_sat_e_upper = nan;
            fit_iv_curve_out.I_sat_e = nan;
            fit_iv_curve_out.I_sat_e_lower = nan;

            fit_iv_curve_out.I_sat_i = nan;
            fit_iv_curve_out.I_sat_i_upper = nan;
            fit_iv_curve_out.I_sat_i_lower = nan;
            
            %Initialize the diagonstics output structure
            diagnostics = struct;
            diagnostics.error_msg = [];
            
            diagnostics.p_I = [];
            diagnostics.p_DI = [];
            diagnostics.p_D2I = [];
            
            diagnostics.logI = [];
            diagnostics.p_logI = [];
            diagnostics.p_DlogI = [];
            
            diagnostics.V_tf = [];
            diagnostics.I_tf = [];
            
            diagnostics.n_e_expected = [];
            diagnostics.n_i_expected = [];
            diagnostics.T_e_expected = [];
            diagnostics.T_i_expected = [];
            
            diagnostics.I_sat_i_expected = [];
            diagnostics.I_sat_e_expected = [];

            %FUNCTION PARAMETERS
            %Frame length fraction for Savitzky-Golay filter
            frame_length_fraction = 1/4
            %Polynomial grade for Savitzky-Golay filter
            n_sgolay = 5
        %     %Plasma detector: max/abs(min) threshold
        %     max_over_abs_min_th = 4;
        %     %Plasma detector: max(abs)/delta_current threshold
        %     max_abs_over_delta_current_th = 4;
            %Multiplication factor for the up-shifting of logI
            min_multiply_factor = 1.1
            %Polynomial grade for the approximation of I fixed to 5
            %Polynomial grade for the approximation of logI fixed to 5
        %     %Maximum Temperature allowed [eV]
        %     T_e_max = 30
        %     %Minimum Temperature allowed [eV]
        %     T_e_min = 1
        %     %Minimum absolute magnitude of Ion Saturation Current allowed [A]
        %     I_sat_i_abs_min = 0.001; %E-GUN PLASMA
        %     %I_sat_i_abs_min = 0.0001; %ECRH PLASMA
            %Voltage where to pick the ion saturation current [V]
            V_i_sat_leftward_offset = 100

            L = length(V_axis);
            frame_length = round(frame_length_fraction*L);
            if mod(frame_length,2) == 0
                frame_length = frame_length + 1;
            end

            %Apply Savitzky–Golay filter
            I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
            fprintf('Savitzky–Golay filter applied\n');

%             %Plasma detector:
%             if (max(I_of_V_sgolay)/abs(min(I_of_V_sgolay)) < max_over_abs_min_th && max(abs(I_of_V_sgolay))/delta_I < max_abs_over_delta_current_th)
%                 fprintf('No plasma detected: abort parameter extraction\n');
%                 %Exit from the current iteration in 1:NC_tot
%                 %and go to the next
%                 continue;
%             else
%                 fprintf('Plasma detected\n');
%             end

            %Polynomial approximations of I(V), I'(V), I''(V)
            n_poly = 5;
            p_I = polyfit(V_axis,I_of_V_sgolay,n_poly);
            p_DI = polyder(p_I);
            p_D2I = polyder(p_DI);

            %Shift the I = I(V) curve upward to avoid negative values in the logarithm
            if min(I_of_V_sgolay) < 0
                shift = -min_multiply_factor*min(I_of_V_sgolay);
                fprintf(['I(V) shifted upward, I_shift = ',num2str(shift),' A\n']);
            else
                shift = 0;
            end

            %Polynomial approximation of log(I+I_shift) and log(I+I_shift)'
            logI = log(I_of_V_sgolay + shift);
            fprintf('Logarithm of I=I(V) + I_shift taken\n');
            n_poly_log = 5;
            p_logI = polyfit(V_axis,logI,n_poly_log);
            fprintf(['Logarithm of I=I(V) + I_shift approximated with polynomial of degree ',num2str(n_poly_log),'\n']);
            p_DlogI = polyder(p_logI);
            fprintf('Polynomial derivative of log(I+I_shfit) taken\n');

%             %Extract floating potential
%             [ ~, ind_V_float ] = min(abs(I_of_V_sgolay));
%             if ismember(ind_V_float,1:length(V_axis))
%                 V_float = V_axis(ind_V_float);
%                 fprintf(['Floating potential extracted: V_float = ',num2str(V_float),' V\n']);
%             else
%                 fprintf('Floating potential is out of range: plasma parameter extraction aborted\n');
%                 continue
%             end

            %Extract 1st approx of plasma potential (a root of D2I)
            r = roots(p_D2I);
            r = r(r > min(V_axis) & r < max(V_axis) & r >= 0);
            if isempty(r)
                fprintf('I''''(V) has no roots greater than V_float: fitIvCurve aborted\n');
            end

            if isempty(r)
                fprintf('Plasma potential is out of range\n');
                return
            end

            V_plasma = min(r);
            [~,V_plasma_ind] = min(abs(V_axis-V_plasma));


            %Look for peaks of DLogI
            DlogI = polyval(p_DlogI,V_axis);
            pks = findpeaks(DlogI);
            if ~isempty(pks)
                mp = pks(1);
            else
                fprintf('Derivative of log(I+I_shift) has no peaks\n');
                mp = 1/10;
            end

            %NEW FIT RANGE: from beginning of V_axis to plasma potential
            fit_range = 1:V_plasma_ind;
            if isempty(fit_range)
                fprintf('Fit range does not contain any indices \n');
            else
                fprintf('Fit range indentified\n');
            end

            V_tf = V_axis(fit_range);
            I_tf = I_of_V_sgolay(fit_range);

            fprintf('Starting fitting procedure...\n');
            [xData, yData] = prepareCurveData( V_tf, I_tf );

            %FIT PARAMETERS (A LOT)
            ft = fittype( 'I_0 + a*V + b*exp(c*V)', 'independent', 'V', 'dependent', 'I' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [-Inf 0 0.001 0.01];
            opts.MaxFunEvals = 600;
            opts.MaxIter = 600;
            opts.Robust = 'Bisquare';
            opts.StartPoint = [0 0.0001 0.1 mp];
            opts.TolFun = 0.0001;
            opts.TolX = 0.0001;
            opts.Upper = [0 Inf 1 Inf];

            % Fit model to data.
            try
                [fitresult, gof] = fit( xData, yData, ft, opts );
                fit_ended_well = true;
            catch ex
                fprintf('Fit procedure failed: fitIvCurve aborted');
                fit_ended_well = false;
            end
            if fit_ended_well
                cv = coeffvalues(fitresult);
                ci = confint(fitresult);
                fprintf('Fitting procedure successful!\n');

                I_0 = cv(1);
                a = cv(2);
                b = cv(3);
                c = cv(4);

                %Fitted model:
                I_of_V_fitted = I_0 + a.*V_axis + b.*exp(c.*V_axis);

                %Plasma parameters:

                %Floating potential
                %Intersect between fitted I(V) curve and horizontal axis
                [ ~, V_float_ind ] = min(abs(I_of_V_fitted));
                if ismember(V_float_ind,1:length(V_axis))
                    V_float = V_axis(V_float_ind)
                else
                    fprintf('Floating potential is out of range\n');
                    V_float = nan;
                end

                %Try to adjust plasma potential value
                %Midpoint between a root of D2I and a root of D3I
                p_D3I = polyder(p_D2I);
                r = roots(p_D3I);

                if isempty(r)
                    fprintf('Plasma potential adjustment not possible.\n');
                else
                    V_plasma = (V_plasma + max(r))/2;
                    [~,V_plasma_ind] = min(abs(V_axis-V_plasma));
                end
                fprintf(['Plasma potential extracted: V_plasma = ',num2str(V_plasma),' V\n']);


                %Electron Temperature
                %Inverse of the exponential coeffcicient
                T_fit = 1/c;
                fprintf(['Electron temperature extracted: T_e = ',num2str(T_fit),' eV\n']);
                T_upper = 1/ci(1,4);
                if T_upper < T_fit
                    T_upper = T_fit;
                end
                T_lower = 1/ci(2,4);
                if T_lower > T_fit
                    T_lower = T_fit;
                end              


                %Ion saturation current
                %Pure ion current somewhere far in the negative bias voltage region
                a_fit = a;
                I_sat_i_fit = I_0 + a_fit*(V_float - V_i_sat_leftward_offset);
                fprintf(['Ion saturation current extracted: I_sat_i_fit = ',num2str(I_sat_i_fit),' A\n']);

                a_upper = ci(1,1);
                if a_upper < a_fit
                    a_upper = a_fit;
                end
                I_sat_i_upper = I_0 + a_upper*(V_float - V_i_sat_leftward_offset);
                a_lower = ci(2,1);
                if a_lower > a_fit
                    a_lower = a_fit;
                end
                I_sat_i_lower = I_0 + a_lower*(V_float - V_i_sat_leftward_offset);

    %             if T_fit > T_e_min && T_fit < T_e_max
    %                 fprintf('Electron temperature is in range\n');
    %             else
    %                 fprintf('Temperature out of range.\n');
    %             end

    %             if abs(I_sat_i_fit) > I_sat_i_abs_min
    %                 fprintf('Ion Saturation Current is in range.\n');
    %             else
    %                 fprintf('Ion Saturation Current is out of range.\n');
    %             end

                %If arrived to this point, it means that parameter extraction
                %went good during this iteration: assign the extracted parameters to the output
                %structure
                fprintf('Parameter extraction successful, assigning parameters values to output structure...\n');
            end
            if fit_ended_well
                    fit_iv_curve_out.I_of_V_sgolay = I_of_V_sgolay;
                    fit_iv_curve_out.fitresult = fitresult;
                    fit_iv_curve_out.gof = gof;
                    fit_iv_curve_out.I_of_V_fitted = I_0 + a.*V_axis + b.*exp(c.*V_axis);
                    fit_iv_curve_out.I_ion_fitted = I_0 + a.*V_axis;

                    fit_iv_curve_out.T_e = T_fit;
                    fit_iv_curve_out.T_e_upper = T_upper;
                    fit_iv_curve_out.T_e_lower = T_lower;

                    fit_iv_curve_out.V_float = V_float;
                    fit_iv_curve_out.V_plasma = V_plasma;

                    fit_iv_curve_out.I_sat_i = I_sat_i_fit;
                    fit_iv_curve_out.I_sat_i_upper = I_sat_i_upper;
                    fit_iv_curve_out.I_sat_i_lower = I_sat_i_lower;
            else
                    fit_iv_curve_out.I_of_V_sgolay = I_of_V_sgolay;
                    fit_iv_curve_out.fitresult_cell_array = [];
                    fit_iv_curve_out.gof_cell_array = [];
                    fit_iv_curve_out.I_of_V_fitted = [];
                    fit_iv_curve_out.I_ion_fitted = [];

                    fit_iv_curve_out.T_e = nan;
                    fit_iv_curve_out.T_e_upper = nan;
                    fit_iv_curve_out.T_e_lower = nan;

                    fit_iv_curve_out.V_float = nan;
                    fit_iv_curve_out.V_plasma = nan;

                    fit_iv_curve_out.I_sat_i = nan;
                    fit_iv_curve_out.I_sat_i_upper = nan;
                    fit_iv_curve_out.I_sat_i_lower = nan;
            end
            fprintf('Done\n');
            %If arrived untill here it means that everything went good:
            %set the return status properly
            fprintf('\n\nfitIvCurveVers01 executed successfully\n');
            return_status = 1;
        end
        
        %% VERSION 2 - STABLE - MUTUAL CONSISTENCY OF PARAMETERS - SIMULTANEOUS FIT OF I_sat_i, I_sat_e, V_plasma, V_float (no sheath expansion effect)
        function [return_status, fit_iv_curve_out, diagnostics] = fitIvCurveVers02(V_axis,I_of_V,S_probe,A_ion,options)
            
            if isfield(options,'diag_msg_flag')
                diag_msg_flag = options.diag_msg_flag;
            else
                diag_msg_flag = false;
            end
            
            %Expected values (based on empirical observation)
            %Electron density
            if isfield(options,'n_e_expected')
                n_e_expected = options.n_e_expected;
            else
                n_e_expected = 1e19;
            end
            %Ion density
            n_i_expected = n_e_expected;
            
            %Electron temperature
            if isfield(options,'T_e_expected')
                T_e_expected = options.T_e_expected;
            else
                T_e_expected = 10;
            end
            
            %Max allowed electron temperature
            if isfield(options,'T_e_max_allowed')
                T_e_max_allowed = options.T_e_max_allowed;
            else
                T_e_max_allowed = 30;
            end

            
            if diag_msg_flag
                fprintf('Running fitIvCurveVers02 with following parameters:\n\n');
            end
            
            if ~iscolumn(V_axis)
                V_axis = V_axis';
            end
            if ~iscolumn(I_of_V)
                I_of_V = I_of_V';
            end
            
            %Prepare return status flag
            return_status = 0;
            %Initialize the fit_iv_curve_out output structure
            fit_iv_curve_out.fitresult_cell_array = [];
            fit_iv_curve_out.gof_cell_array = [];
            fit_iv_curve_out.I_of_V_fitted = [];
            fit_iv_curve_out.I_ion_fitted = [];

            fit_iv_curve_out.T_e = nan;
            fit_iv_curve_out.T_e_upper = nan;
            fit_iv_curve_out.T_e_lower = nan;

            fit_iv_curve_out.V_float = nan;
            fit_iv_curve_out.V_plasma = nan;

            fit_iv_curve_out.I_sat_i = nan;
            fit_iv_curve_out.I_sat_i_upper = nan;
            fit_iv_curve_out.I_sat_i_lower = nan;
            
            %Initialize the diagonstics output structure
            diagnostics = struct;
            diagnostics.error_msg = [];
            
            diagnostics.p_I = [];
            diagnostics.p_DI = [];
            diagnostics.p_D2I = [];
            
            diagnostics.logI = [];
            diagnostics.p_logI = [];
            diagnostics.p_DlogI = [];
            
            diagnostics.V_tf = [];
            diagnostics.I_tf = [];
            
            diagnostics.n_e_expected = [];
            diagnostics.n_i_expected = [];
            diagnostics.T_e_expected = [];
            diagnostics.T_i_expected = [];
            
            diagnostics.I_sat_i_expected = [];
            diagnostics.I_sat_e_expected = [];
            
            %CONVENTION:
            %ELECTRON CURRENT IS POSITIVE
            %ION CURRENT IS NEGATIVE

            %FUNCTION PARAMETERS
            %Frame length fraction for Savitzky-Golay filter
            frame_length_fraction = 1/4;
            %Polynomial grade for Savitzky-Golay filter
            n_sgolay = 5;
            %Multiplication factor for the up-shifting of logI
            min_multiply_factor = 1.1;
            %Polynomial grade for the approximation of I fixed to 5
            %Polynomial grade for the approximation of logI fixed to 5
            L = length(V_axis);
            frame_length = round(frame_length_fraction*L);
            if mod(frame_length,2) == 0
                frame_length = frame_length + 1;
            end
            %Apply Savitzky–Golay filter
            I_of_V_sgolay = sgolayfilt(I_of_V,n_sgolay,frame_length);
            diagnostics.I_of_V_sgolay = I_of_V_sgolay;
            if diag_msg_flag
                fprintf('Savitzky–Golay filter applied\n'); %#ok<UNRCH>
            end
            
            %Spline fit of I(V)
            spline_I = fit(V_axis,I_of_V,'smoothingspline','SmoothingParam',1e-4);
            I_of_V_spline = ppval(spline_I.p,V_axis);
            diagnostics.I_of_V_spline = I_of_V_spline;
            
            %Piecewise polynomial approximations of I(V), I'(V), I''(V)
            diagnostics.pp_I = spline_I.p;
           
            %Polynomial approximations of I(V), I'(V), I''(V)
            n_poly = 5;
            ws = warning('off','all');
            p_I = polyfit(V_axis,I_of_V_sgolay,n_poly);
            warning(ws)
            p_DI = polyder(p_I);
            p_D2I = polyder(p_DI);
            
            diagnostics.p_I = p_I;
            diagnostics.p_DI = p_DI;
            diagnostics.p_D2I = p_D2I;

            %Shift the I = I(V) curve upward to avoid negative values in the logarithm
            if min(I_of_V_sgolay) < 0
                shift = -min_multiply_factor*min(I_of_V_sgolay);
                if diag_msg_flag
                    fprintf(['I(V) shifted upward, I_shift = ',num2str(shift),' A\n']); %#ok<UNRCH>
                end
            else
                shift = 0;
            end

            %Polynomial approximation of log(I+I_shift) and log(I+I_shift)'
            logI = log(I_of_V_sgolay + shift);
            if diag_msg_flag
                fprintf('Logarithm of I=I(V) + I_shift taken\n');
            end
            n_poly_log = 5;
            ws = warning('off','all');
            p_logI = polyfit(V_axis,logI,n_poly_log);
            warning(ws)
            if diag_msg_flag
                fprintf(['Logarithm of I=I(V) + I_shift approximated with polynomial of degree ',num2str(n_poly_log),'\n']);
            end
            p_DlogI = polyder(p_logI);
            if diag_msg_flag
                fprintf('Polynomial derivative of log(I+I_shfit) taken\n');
            end
            
            diagnostics.logI = logI;
            diagnostics.p_logI = p_logI;
            diagnostics.p_DlogI = p_DlogI;
            
            
            

            %Plasma potential (a root of D2I) (center of exponential region)
            r_arr = sort(roots(p_D2I));
            V_plasma = [];
            try
                V_plasma = r_arr(2);
            catch
                if diag_msg_flag
                error_msg = 'Plasma potential could not be found as root of D2I';
                diagonstics.error_msg = error_msg;
                fprintf(strcat(diagnostics.error_msg,'\n'));
                end

                return
            end

            [~,V_plasma_ind] = min(abs(V_axis-V_plasma));


            %Look for peaks of DLogI
            DlogI = polyval(p_DlogI,V_axis);
            pks = findpeaks(DlogI);
            if ~isempty(pks)
                inv_Temp_1st_approx = pks(1);
            else
                if diag_msg_flag
                    error_msg = 'Derivative of log(I+I_shift) has no peaks: fitIvCurve aborted';
                    diagnostics.error_msg = error_msg;
                    fprintf(strcat(diagnostics.error_msg,'\n'));
                end

                return
            end

            %FIT RANGE: from beginning of V_axis to one step right to plasma potential
            fit_range = 1:(V_plasma_ind+1);
            if isempty(fit_range) || fit_range(end) > length(V_axis)
                if diag_msg_flag
                    error_msg = 'Fit range is invalid';
                    diagonstics.error_msg = error_msg;
                    fprintf(strcat(diagnostics.error_msg,'\n'));
                end
                
                return
            else
                if diag_msg_flag
                    fprintf('Fit range indentified\n');
                end
            end
            

            V_tf = V_axis(fit_range);
            %I_tf = I_of_V_sgolay(fit_range);
            %I_tf = I_of_V_spline(fit_range);
            I_tf = I_of_V(fit_range);
            diagnostics.V_tf = V_tf;
            diagnostics.I_tf = I_tf;
            

            
            T_i_expected = 300*(PlasmaPhysics.k_B/abs(PlasmaPhysics.e));
            diagnostics.n_e_expected = n_e_expected;
            diagnostics.n_i_expected = n_i_expected;
            diagnostics.T_e_expected = T_e_expected;
            diagnostics.T_i_expected = T_i_expected;
            
            M_i = PlasmaPhysics.AMU*A_ion;

            I_sat_i_expected = (1/4)*PlasmaPhysics.e*n_i_expected*(1.6)*sqrt(abs(PlasmaPhysics.e)*T_i_expected/M_i)*S_probe;
            %I_sat_e_expected = -(1/4)*PlasmaPhysics.e*n_e_expected*(1.6)*sqrt(abs(PlasmaPhysics.e)*T_e_expected/PlasmaPhysics.m_e)*S_probe;
            I_sat_e_expected = I_of_V(V_plasma_ind+1);
            
            diagnostics.I_sat_i_expected = I_sat_i_expected;
            diagnostics.I_sat_e_expected = I_sat_e_expected;

            if diag_msg_flag
                fprintf('Starting fitting procedure...\n');
            end
            
            [xData, yData] = prepareCurveData( V_tf, I_tf );
            
            inv_T_e_expected = inv_Temp_1st_approx;
            
            % Set up fittype and options.
            ft = fittype( strcat('I_sat_i + I_sat_e*exp((V-V_plasma)*inv_T_e)'), 'independent', 'V', 'dependent', 'I' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            
            opts.StartPoint = [I_sat_e_expected I_sat_i_expected V_plasma inv_T_e_expected];
            
            opts.Upper = [I_sat_e_expected       0    V_plasma+100    Inf ];
            opts.Lower = [I_sat_e_expected      -1    V_plasma-100    1/T_e_max_allowed];

            % Fit model to data.
            try
                [fitresult, gof] = fit( xData, yData, ft, opts);
                fit_ended_well = true;
            catch err
                if diag_msg_flag
                    disp(err.message);
                    fprintf('Fit procedure failed: fitIvCurve aborted');
                end
                fit_ended_well = false;
            end
            
            if fit_ended_well
                cv = coeffvalues(fitresult);
                ci = confint(fitresult);
                if diag_msg_flag
                    fprintf('Fitting procedure successful!\n');
                    disp(fitresult);
                    fprintf('Goodness:\n');
                    disp(gof);
                end

                I_sat_e_upper = ci(2,1);
                I_sat_e = cv(1);
                I_sat_e_lower = ci(1,1);

                I_sat_i_upper = ci(2,2);
                I_sat_i = cv(2);
                I_sat_i_lower = ci(1,2);
                
                V_plasma_upper = ci(2,3);
                V_plasma = cv(3);
                V_plasma_lower = ci(1,3);
                
                inv_T_e_upper = ci(2,4);
                inv_T_e = cv(4);
                inv_T_e_lower = ci(1,4);
                
                T_e = 1/inv_T_e;

                %Fitted model:
                I_of_V_fitted = I_sat_i + I_sat_e.*exp((V_axis-V_plasma).*inv_T_e);
                
                %Floating potential
                %Intersect between fitted I(V) curve and horizontal axis
                delta_V = min(abs(diff(V_axis)));
                V_axis_finer = min(V_axis):delta_V:max(V_axis);
                I_of_V_finer = feval(fitresult,V_axis_finer);
                [~,V_float_ind] = min(abs(I_of_V_finer));
                V_float = V_axis_finer(V_float_ind);


                %Plasma parameters:
                if diag_msg_flag
                fprintf(['Floating potential: V_float = ',num2str(V_float),' V\n']);
                fprintf(['Plasma potential: V_plasma = ',num2str(V_plasma),' V\n']);
                fprintf('Confidence bounds:\n');
                fprintf(['V_plasma_upper = ',num2str(V_plasma_upper),' V\n']);
                fprintf(['V_plasma_lower = ',num2str(V_plasma_lower),' V\n']);
                fprintf(['Electron temperature: T_e = ',num2str(T_e),' eV\n']);
                fprintf('Confidence bounds:\n');
                fprintf(['T_e_upper = ',num2str(1/inv_T_e_lower),' eV\n']);
                fprintf(['T_e_lower = ',num2str(1/inv_T_e_upper),' eV\n']);
                %fprintf(['Electron saturation current: I_sat_e = ',num2str(I_sat_e),' A\n']);
                fprintf(['Ion saturation current: I_sat_i = ',num2str(I_sat_i),' A\n']);
                fprintf('Confidence bounds:\n');
                fprintf(['I_sat_i_upper = ',num2str(I_sat_i_upper),' A\n']);
                fprintf(['I_sat_i_lower = ',num2str(I_sat_i_lower),' A\n']);
                end

                %If arrived to this point, it means that parameter extraction
                %went good during this iteration: assign the extracted parameters to the output
                %structure
                if diag_msg_flag
                    fprintf('Parameter extraction successful, assigning parameters values to output structure...\n');
                end

                fit_iv_curve_out.fitresult = fitresult;
                fit_iv_curve_out.gof = gof;
                fit_iv_curve_out.I_of_V_fitted = I_of_V_fitted;

                fit_iv_curve_out.V_float = V_float;
                
                fit_iv_curve_out.V_plasma_upper = V_plasma_upper;
                fit_iv_curve_out.V_plasma = V_plasma;
                fit_iv_curve_out.V_plasma_lower = V_plasma_lower;
                
                fit_iv_curve_out.I_sat_e_upper = I_sat_e_upper;
                fit_iv_curve_out.I_sat_e = I_sat_e;
                fit_iv_curve_out.I_sat_e_lower = I_sat_e_lower;

                fit_iv_curve_out.I_sat_i_upper = I_sat_i_upper;
                fit_iv_curve_out.I_sat_i = I_sat_i;
                fit_iv_curve_out.I_sat_i_lower = I_sat_i_lower;
                
                fit_iv_curve_out.T_e_upper = 1/inv_T_e_lower;
                fit_iv_curve_out.T_e = 1/inv_T_e;
                fit_iv_curve_out.T_e_lower = 1/inv_T_e_upper;
                
                if diag_msg_flag
                    fprintf('Done\n');
                    %If arrived untill here it means that everything went good:
                    %set the return status properly
                    fprintf('\n\nfitIvCurveVers02 executed successfully\n');
                end
                return_status = 1;
            else
                
                if diag_msg_flag
                    fprintf('Done\n');
                    %If arrived untill here it means that everything went good:
                    %set the return status properly
                    fprintf('\n\nfitIvCurveVers02 executed, but no parameters extracted\n');
                end
                return_status = 0;
            end
            
        end
        
        %% VERSION 3 - STABLE - Sheath expansion effect and precise determination of V_plasma, I_sat_e, I_sat_i, T_e
        function [return_status, out, diagnostics] = fitIvCurveVers03(V_axis,I_of_V,S_probe,A_ion,options)
            
            %Prepare output structures
            out = struct;
            diagnostics = struct;
            
            if isfield(options,'V_bars')
                V_bars = options.V_bars;
                if ~iscolumn(V_bars)
                    V_bars = V_bars';
                end
            else
                V_bars = zeros(length(V_axis),1);
            end
            out.V_bars = V_bars;
            
            if isfield(options,'I_bars')
                I_bars = options.I_bars;
                if ~iscolumn(I_bars)
                    I_bars = I_bars';
                end
            else
                I_bars = zeros(length(V_axis),1);
            end
            out.I_bars = I_bars;
            
            if isfield(options,'V_offset_correction')
                V_axis = V_axis + options.V_offset_correction;
            end
            if isfield(options,'I_offset_correction')
                I_of_V = I_of_V + options.I_offset_correction;
            end
            
            
            if isfield(options,'diag_msg_flag')
                diag_msg_flag = options.diag_msg_flag;
            else
                diag_msg_flag = false;
            end
            
            %Expected values (based on empirical observation)
            %Electron density
            if isfield(options,'n_e_expected')
                n_e_expected = options.n_e_expected;
            else
                n_e_expected = 1e19;
            end
            %Ion density
            n_i_expected = n_e_expected;
            
            %Electron temperature
            if isfield(options,'T_e_expected')
                T_e_expected = options.T_e_expected;
            else
                T_e_expected = 10;
            end
            
            %Max allowed electron temperature
            if isfield(options,'T_e_max_allowed')
                T_e_max_allowed = options.T_e_max_allowed;
            else
                T_e_max_allowed = 30;
            end

            
            if diag_msg_flag
                fprintf('Running fitIvCurveVers03\n\n');
            end
            
            if ~iscolumn(V_axis)
                V_axis = V_axis';
            end
            if ~iscolumn(I_of_V)
                I_of_V = I_of_V';
            end
            
            %Prepare return status flag
            return_status = 0;
            

            
            out.V_axis = V_axis;
            out.I_of_V = I_of_V;
            
            %CONVENTION:
            %ELECTRON CURRENT IS POSITIVE
            %ION CURRENT IS NEGATIVE
            
            % 0) Approximate dataset with a spline
            I_of_V_pp = spline(V_axis,I_of_V);
            out.I_of_V_pp = I_of_V_pp;
            %Shift Up and left
            I_of_V_up = I_of_V+I_bars;
            I_of_V_up_pp = spline(V_axis,I_of_V_up);
            out.I_of_V_up_pp = I_of_V_up_pp;
            r_upleft = fnzeros(I_of_V_up_pp);
            %Shift Down and right
            I_of_V_down = I_of_V-I_bars;
            I_of_V_down_pp = spline(V_axis,I_of_V_down);
            out.I_of_V_down_pp = I_of_V_down_pp;
            r_downright = fnzeros(I_of_V_down_pp);
            
            % 1) Find I = 0 Voltage and V_float (ideally they should
            % coincide, they will NOT coincide if hysteresis is present)
            
            [~,V_null_current] = min(abs(I_of_V));
            r = fnzeros(I_of_V_pp);
            V_null_current = r(1);             
            V_null_current_bounds(1) = r_upleft(1)-V_bars(1);
            V_null_current_bounds(2) = r_downright(1)+V_bars(1);
            out.V_null_current = V_null_current;
            out.V_null_current_bounds = V_null_current_bounds;
            
            V_float_estimate_bounds = [nan,nan];
            % If it is provided as input parameter...
            if isfield(options,'include_external_V_float') && options.include_external_V_float
                % Take the average
                V_float_estimate = mean(options.external_V_float);
                [~,V_float_estimate_ind] = min(abs(V_axis-V_float_estimate));
                V_float_std = std(options.external_V_float);
                V_float_estimate_bounds(1) = V_float_estimate - V_float_std;
                V_float_estimate_bounds(2) = V_float_estimate + V_float_std;
            else
                %Otherwise...
                %...(intercept of I_of_V with I = 0 line)
                [~,V_float_estimate_ind] = min(abs(I_of_V));
                V_float_estimate = V_null_current;          
                V_float_estimate_bounds(1) = V_null_current_bounds(1);
                V_float_estimate_bounds(2) = V_null_current_bounds(2);
            end
            out.V_float_estimate = V_float_estimate;
            out.V_float_estimate_ind = V_float_estimate_ind;
            out.V_float_estimate_bounds = V_float_estimate_bounds;
            
            % 2) Linear fit of ion current due to sheath expansion:
            %   I_i(V) = I_i_0 + a*V in the interval [min(V_axis),(almost)V_float] 
            
            fit_range_I_i = 1:V_float_estimate_ind;
            
            V_tf_I_i = V_axis(fit_range_I_i);
            I_tf_I_i = I_of_V(fit_range_I_i);
            
            [V_tf_exp_Data, I_e_tf_exp] = prepareCurveData( V_tf_I_i, I_tf_I_i );

            % Set up fittype and options.
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            opts.Normalize = 'off';
            opts.Robust = 'Bisquare';

            % Fit model to data.
            [fitresult_I_i, gof_I_i] = fit( V_tf_exp_Data, I_e_tf_exp, ft, opts );
            out.fitresult_I_i = fitresult_I_i;
            out.gof_I_i = gof_I_i;
            ci= confint(fitresult_I_i);
            
            I_i_offset = fitresult_I_i.p2;
            I_i_offset_bounds = ci(:,2);
            out.I_i_offset = I_i_offset;
            I_i_offset_bounds = sort(I_i_offset_bounds);
            I_i_offset_bounds(1) = I_i_offset_bounds(1) - I_bars(1);
            I_i_offset_bounds(2) = I_i_offset_bounds(2) + I_bars(1);
            out.I_i_offset_bounds = I_i_offset_bounds;
            
            sheath_coeff = fitresult_I_i.p1;
            sheath_coeff_bounds = ci(:,1);
            out.sheath_coeff = sheath_coeff;
            sheath_coeff_bounds = sort(sheath_coeff_bounds);
            out.sheath_coeff_bounds = sheath_coeff_bounds;
            
            I_i = fitresult_I_i.p2 + sheath_coeff.*V_axis;
            out.I_i = I_i;
            
            % 3) Focus on electron current I_e = I_of_V - I_i
            I_e = I_of_V - I_i;
            out.I_e = I_e;
            

            % 4) Logarithm of I_e and localization of the knee
            log_I_e = log(I_e);
            out.log_I_e = log_I_e;
            
            positive_range = (V_float_estimate_ind+1):length(V_axis);
            positive_range_D = positive_range(1:end-1);
            positive_range_D2 = positive_range_D(1:end-1);
            out.positive_range = positive_range;
            out.positive_range_D = positive_range_D;
            out.positive_range_D2 = positive_range_D2;
            
            D_log_I_e = diff(log_I_e);
            D_log_I_e_smooth = smoothdata(D_log_I_e,'rlowess',10);
            D2_log_I_e = diff(D_log_I_e_smooth);
            D2_log_I_e_smooth = smoothdata(D2_log_I_e,'rlowess',10);
            out.D_log_I_e = D_log_I_e;
            out.D_log_I_e_smooth = D_log_I_e_smooth;
            out.D2_log_I_e = D2_log_I_e;
            out.D2_log_I_e_smooth = D2_log_I_e_smooth;
            
            %[~,knee_ind_in_positive_range] = min(D2_log_I_e_smooth(positive_range_D2));
            [~,knee_ind_in_positive_range] = min(D2_log_I_e(positive_range_D2));
            knee_ind = positive_range(1) + knee_ind_in_positive_range;
            out.knee_ind_in_positive_range = knee_ind_in_positive_range;
            out.knee_ind = knee_ind;
            
            V_inflection = V_axis(knee_ind);
            out.V_inflection = V_inflection;
            
            % 5) Slope of linear fit of log(I_e) before the knee gives the inverse of T_e
            V_tf = V_axis(positive_range(1):knee_ind);
            log_I_e_tf = log_I_e(positive_range(1):knee_ind);
            %Linear region in semilog
            [V_data, log_I_e_data] = prepareCurveData( V_tf, log_I_e_tf );
            % Set up fittype and options.
            ft = fittype( 'poly1' );
            opts = fitoptions( 'Method', 'LinearLeastSquares' );
            opts.Robust = 'Bisquare';
            % Fit model to data.
            [fitresult_log_I_e, gof_log_I_e] = fit( V_data, log_I_e_data, ft, opts );
            out.fitresult_log_I_e = fitresult_log_I_e;
            out.gof_log_I_e = gof_log_I_e;
            ci = confint(fitresult_log_I_e);
            
            inv_T_e_estimate = fitresult_log_I_e.p1;
            inv_T_e_estimate_bounds = ci(:,1);
            T_e_estimate = 1/inv_T_e_estimate;
            T_e_estimate_bounds = (1./inv_T_e_estimate_bounds)';
            T_e_estimate_bounds = sort(T_e_estimate_bounds); %#ok<TRSRT>
            
            out.inv_T_e_estimate = inv_T_e_estimate;
            out.inv_T_e_estimate_bounds = inv_T_e_estimate_bounds;
            out.T_e_estimate = T_e_estimate;
            out.T_e_estimate_bounds = T_e_estimate_bounds;
            
            % 6) V_plasma is calculated from V_float_estimate and T_e
            m_e = PlasmaPhysics.m_e;
            m_i = A_ion*PlasmaPhysics.AMU;
            out.m_e = m_e;
            out.m_i = m_i;
            
            V_plasma_estimate = V_float_estimate - T_e_estimate*log(0.6*sqrt(2*pi*m_e/m_i));
            out.V_plasma_estimate = V_plasma_estimate;
            
            V_plasma_estimate_bounds = [nan,nan];
            V_plasma_estimate_bounds(1) = min(V_float_estimate_bounds) - min(T_e_estimate_bounds)*log(0.6*sqrt(2*pi*m_e/m_i));
            V_plasma_estimate_bounds(2) = max(V_float_estimate_bounds) - max(T_e_estimate_bounds)*log(0.6*sqrt(2*pi*m_e/m_i));
            V_plasma_estimate_bounds = sort(V_plasma_estimate_bounds);
            out.V_plasma_estimate_bounds = V_plasma_estimate_bounds;
            
            [~,V_plasma_estimate_ind] = min(abs(V_axis - V_plasma_estimate));
            out.V_plasma_estimate_ind = V_plasma_estimate_ind;
            
            V_plasma_estimate_ind_bounds = [nan,nan];
            [~,V_plasma_estimate_ind_bounds(1)] = min(abs(V_axis - V_plasma_estimate_bounds(1)));
            [~,V_plasma_estimate_ind_bounds(2)] = min(abs(V_axis - V_plasma_estimate_bounds(2)));
            V_plasma_estimate_ind_bounds = sort(V_plasma_estimate_ind_bounds);
            out.V_plasma_estimate_ind_bounds = V_plasma_estimate_ind_bounds;
                
%             V_plasma_estimate = V_inflection;
%             out.V_plasma_estimate = V_plasma_estimate;
%             V_plasma_estimate_bounds = [V_plasma_estimate - V_bars(1),V_plasma_estimate + V_bars(1)];
%             out.V_plasma_estimate_bounds = V_plasma_estimate_bounds;
%             
%             V_plasma_estimate_ind_bounds = [nan,nan];
%             [~,V_plasma_estimate_ind_bounds(1)] = min(abs(V_axis - V_plasma_estimate_bounds(1)));
%             [~,V_plasma_estimate_ind_bounds(2)] = min(abs(V_axis - V_plasma_estimate_bounds(2)));
%             V_plasma_estimate_ind_bounds = sort(V_plasma_estimate_ind_bounds);
%             out.V_plasma_estimate_ind_bounds = V_plasma_estimate_ind_bounds;
            
%             %Recalculate V_float
%             V_float_estimate = V_plasma_estimate + T_e_estimate*log(0.6*sqrt(2*pi*m_e/m_i));
%             out.V_float_estimate = V_float_estimate;
%             V_float_estimate_bounds = V_plasma_estimate_bounds + T_e_estimate_bounds.*log(0.6*sqrt(2*pi*m_e/m_i));
%             out.V_float_estimate_bounds = V_float_estimate_bounds;
            
            % 7) I_sat_e = I_e(V_plasma) and I_sat_i = I_i(V_plasma) 
            
            I_sat_e_estimate = ppval(I_of_V_pp,V_plasma_estimate);
            out.I_sat_e_estimate = I_sat_e_estimate;
            
            I_sat_e_estimate_bounds = ppval(I_of_V_pp,V_plasma_estimate_bounds);
            out.I_sat_e_estimate_bounds = I_sat_e_estimate_bounds;
            
            I_sat_i_estimate = -(0.6)*sqrt(2*pi)*sqrt(m_e/m_i)*I_sat_e_estimate;
            out.I_sat_i_estimate = I_sat_i_estimate;
            
            I_sat_i_estimate_bounds = -(0.6)*sqrt(2*pi)*sqrt(m_e/m_i).*I_sat_e_estimate_bounds;
            out.I_sat_i_estimate_bounds = I_sat_i_estimate_bounds;
            
            n_e_estimate = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,I_sat_i_estimate,T_e_estimate);
            out.n_e_estimate = n_e_estimate;
            
            n_e_estimate_bounds = [nan,nan];
            %Upper -> max I_sat_i, min T_e
            n_e_estimate_bounds(1) = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,max(I_sat_i_estimate_bounds),min(T_e_estimate_bounds));
            %Lower -> min I_sat_i, max T_e
            n_e_estimate_bounds(2) = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,min(I_sat_i_estimate_bounds),max(T_e_estimate_bounds));
            n_e_estimate_bounds = sort(n_e_estimate_bounds);
            out.n_e_estimate_bounds = n_e_estimate_bounds;
            
%             I_sat_i_estimate = I_i_offset + sheath_coeff.*V_plasma_estimate;
%             out.I_sat_i_estimate = I_sat_i_estimate;
%             
%             I_sat_i_estimate_bounds = [nan,nan];
%             %Min: = maximum offset maximum sheath_coeff maximum V_plasma
%             I_sat_i_estimate_bounds(1) = max(I_i_offset_bounds) + max(sheath_coeff_bounds)*max(V_plasma_estimate_bounds);
%             %Max: = minium offset, minimum sheath_coeff minimum V_plasma
%             I_sat_i_estimate_bounds(2) = min(I_i_offset_bounds) + min(sheath_coeff_bounds)*min(V_plasma_estimate_bounds); 
%             if(max(I_sat_i_estimate_bounds) > 0)
%                 [~,ind] = max(I_sat_i_estimate_bounds);
%                 I_sat_i_estimate_bounds(ind) = 0;
%             end
%             out.I_sat_i_estimate_bounds = I_sat_i_estimate_bounds;
            
            % 8) Exponential fit of I_e in the interval [min(V_axis),V_axis(knee_ind)] giving NEW T_e and I_sat_e
            % fit_range_I_e = 1:knee_ind;
            fit_range_I_e = (V_float_estimate_ind-(knee_ind-V_float_estimate_ind)):knee_ind;
            V_tf_exp = V_axis(fit_range_I_e);
            I_e_tf_exp = I_e(fit_range_I_e);
            out.V_tf_exp = V_tf_exp;
            out.I_e_tf_exp = I_e_tf_exp;
            
            
            [V_tf_exp_Data, I_e_tf_exp_Data] = prepareCurveData( V_tf_exp, I_e_tf_exp );

            % Set up fittype and options.
            ft = fittype( 'exp1' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 0];
            opts.MaxFunEvals = 1000;
            opts.MaxIter = 1000;
            opts.Robust = 'Bisquare';
            opts.StartPoint = [I_sat_e_estimate*exp(-inv_T_e_estimate*V_plasma_estimate) inv_T_e_estimate];

            % Fit model to data.
            [fitresult_I_e, gof_I_e] = fit( V_tf_exp_Data, I_e_tf_exp_Data, ft, opts );
            out.fitresult_I_e = fitresult_I_e;
            out.gof_I_e = gof_I_e;
            
            T_e_fit = 1/(fitresult_I_e.b);
            out.T_e_fit = T_e_fit;
            
            % 9) Rebuild I_fit = I_i_fit + I_e_fit as an abstract function
            I_fit = @(V) fitresult_I_i.p2 + fitresult_I_i.p1.*V + fitresult_I_e.a*exp(V/T_e_fit);
            out.I_fit = I_fit;
            
            %10) Find NEW V_float as the intercept between I_fit and I = 0
            [V_float_fit,fval,exitflag,output] = fzero(I_fit,V_float_estimate);
            out.V_float_fit = V_float_fit;
            
            %11) Find NEW V_plasma as a function of NEW V_float, NEW T_e
            V_plasma_fit = V_float_fit - T_e_fit*log(0.6*sqrt(2*pi*m_e/m_i));
            out.V_plasma_fit = V_plasma_fit;

            %12) Find I_sat_i_fit as I_i(V_plasma)
            I_i = @(V) fitresult_I_i.p2 + fitresult_I_i.p1.*V;
            I_sat_i_fit = I_i(V_plasma_fit);
            out.I_sat_i_fit = I_sat_i_fit;
            
            %12) Find I_sat_e_fit as a funcion of I_sat_i_fit and m_i
            I_sat_e_fit = (0.65)*sqrt(m_i/m_e)*(-I_sat_i_fit);
            out.I_sat_e_fit = I_sat_e_fit;
            
            return_status = 1;

            
        end
        
    end
end

