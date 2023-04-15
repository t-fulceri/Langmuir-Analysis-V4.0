
% addpath('..\..\..\General_Purpose');
% addpath('..\..\..\General_Purpose\periodical_cubic_interpolation');
% addpath('..\..\..\Plasma_Physics');
% addpath('..\..\..\Langmuir Analysis V4.0\Classes');

%Langmuir probe area and gas used
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_probe = L_tip*2*pi*r_tip;
% Ion mass number
A_ion = 39.948;%Argon

keySet = {'stationary','timeseries'};
valueSet = [1 2];
M = containers.Map(keySet,valueSet);

%select = 'timeseries'
select = 'stationary'

%Select rise/fall
rf = 2;

close all;

switch M(select)
    case 1
%         vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_03_07_Langmuir_Vacuum\f_01KHz\';
%         vacuum_filename = 'meas_0002.h5';
%         vacuum_full_file_path = strcat(vacuum_source_folder,vacuum_filename);
%         plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_03_07_Langmuir_ECRH_Central_position_triggered\f_01KHz\';
%         plasma_filename = 'meas_0002.h5';
%         plasma_full_file_path = strcat(plasma_source_folder,plasma_filename);
        
        %plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_05_29_Langmuir_ECRH_Power_Cycle_Characterization\f_500hz_plasma\';
        %plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_06_05_Langmuir_ECRH_Power_Cycle_Characterization_New_Meas_Config\f_500hz_plasma\';
        %plasma_filename = 'meas_0001.h5';
        
        iv_signals_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';
        vf_signal_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_Floating_Potential_Profile\';
        plasma_filename = 'meas_0274.h5';
        
        iv_signals_full_file_path = strcat(iv_signals_source_folder,plasma_filename);
        vf_signal_full_file_path = strcat(vf_signal_source_folder,plasma_filename);
        pp_stationary = PlasmaParamStationary(iv_signals_full_file_path,vf_signal_full_file_path);
        

        switch rf
            case 1
                V_axis = pp_stationary.V_axis_rise;
                I_of_V = pp_stationary.I_of_V_rise;
            case 2
                V_axis = pp_stationary.V_axis_fall;
                I_of_V = pp_stationary.I_of_V_fall;
        end
        
        
        version_selector = 3;
        
        options.include_external_V_float = true;
        
        options.diag_msg_flag = true;
        options.T_e_expected = 10;
        options.T_e_max_allowed = 100;
        options.n_e_expected = 1e17;
        
        pp_stationary.extractPlasmaParameters(S_probe,A_ion,version_selector,options);
        
        output = pp_stationary.fit_iv_curve_out_ca{rf}; %#ok<NASGU>
        save('temp.mat','-struct','output');
        clear output;
        load('temp.mat');
        

        time_axis = pp_stationary.ivt_source.ivs_source_plasma.time_axis;
        V_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.V_signal;
        I_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.I_signal;

        V_signal_vacuum = pp_stationary.ivt_source.ivs_source_vacuum.V_signal;
        I_signal_vacuum = pp_stationary.ivt_source.ivs_source_vacuum.I_signal;

%         V_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.V_signal;
%         I_signal_plasma_CM_corrected = pp_stationary.ivt_source.ivs_source_plasma.I_CM_corrected;
%         I_smooth_exp_fit = pp_stationary.ivt_source.ivs_source_plasma.I_smooth_exp_fit;





        

        
        figure
        subplot(1,3,1);
        title('Semilog plot of electron current');
        xlabel('V [V]');
        ylabel('log(I_e) [a.u.]')
        hold on
        plot(V_axis(positive_range),log_I_e(positive_range),'o');
        plot(V_axis(knee_ind).*[1 1], ylim, '--k','Linewidth',1);
        plot(V_axis(positive_range), fitresult_log_I_e.p2 + inv_T_e_estimate.*V_axis(positive_range), '--r','Linewidth',1);
        ylim([min(log_I_e(positive_range)) max(log_I_e(positive_range))]);
        T_e_string = ['T_{e,est.} = ',num2str(T_e_estimate,2),' eV ','(',num2str(min(T_e_estimate_bounds),2),' eV',',',num2str(max(T_e_estimate_bounds),2),' eV',')'];
        text(V_axis(positive_range(1) + floor(length(positive_range)/4)),log_I_e(positive_range(5)),T_e_string,'FontSize',14);
        hold off
        subplot(1,3,2);
        title('Derivative of logarithm of electron current and its smoothed version');
        xlabel('V [V]');
        ylabel('d(log(I_e))/dV [a.u.]')
        hold on
        plot(V_axis(positive_range_D),D_log_I_e(positive_range_D),'o');
        plot(V_axis(positive_range_D),D_log_I_e_smooth(positive_range_D),'-');
        plot(V_axis(knee_ind).*[1 1], ylim, '--k','Linewidth',1);
        hold off
        subplot(1,3,3);
        title('Second derivative of logarithm of electron current and its smoothed version');
        xlabel('V [V]');
        ylabel('d^2(log(I_e))/dV^2 [a.u.]')
        hold on
        plot(V_axis(positive_range_D2),D2_log_I_e(positive_range_D2),'o');
        plot(V_axis(positive_range_D2),D2_log_I_e_smooth(positive_range_D2),'-');
        plot(V_axis(knee_ind).*[1 1], ylim, '--k','Linewidth',1);
        hold off

        
        figure
        title('First estimation of plasma parameters');
        xlabel('Bias Voltage [V]');
        ylabel('Probe Current [A]');
        hold on
        legend_entries = {};
        

        %errorbar(V_axis,I_of_V,I_bars,'ko');
        errorbar(V_axis,I_of_V,I_bars,I_bars,V_bars,V_bars,'ko')
        legend_entries{end+1} = 'Data points';
        
        plot(V_axis,ppval(I_of_V_pp,V_axis),'-');
        legend_entries{end+1} = 'Spline interpolation';
        
        plot(V_axis,I_i,'-');
        legend_entries{end+1} = 'Ion current';
        
        plot(V_axis,max(I_i_offset_bounds) + max(sheath_coeff_bounds).*V_axis,'--');
        legend_entries{end+1} = 'Ion current Upper';
        
        plot(V_axis,min(I_i_offset_bounds) + min(sheath_coeff_bounds).*V_axis,'--');
        legend_entries{end+1} = 'Ion current Lower';
        
        plot(xlim, [0 0], '-k','Linewidth',1);
        legend_entries{end+1} = 'Zero current line';
        
        plot(V_float_estimate.*[1 1], ylim, '-g','Linewidth',1);
        legend_entries{end+1} = ['V_{float,est.} = ',num2str(V_float_estimate),' V'];
        
        plot(max(V_float_estimate_bounds).*[1 1], ylim, '--g','Linewidth',1);
        legend_entries{end+1} = ['V_{float,est.Upper} = ',num2str(max(V_float_estimate_bounds)),' V'];
        
        plot(min(V_float_estimate_bounds).*[1 1], ylim, '--g','Linewidth',1);
        legend_entries{end+1} = ['V_{float,est.Lower} = ',num2str(min(V_float_estimate_bounds)),' V'];
        
        plot(V_axis(knee_ind).*[1 1], ylim, '-k','Linewidth',1);
        legend_entries{end+1} = ['V_{inflection} = ',num2str(V_axis(knee_ind)),' V'];
        
        plot(V_plasma_estimate.*[1 1], ylim, '-m','Linewidth',1);
        legend_entries{end+1} = ['V_{plasma,est.} = ',num2str(V_plasma_estimate),' V'];
        
        plot(max(V_plasma_estimate_bounds).*[1 1], ylim, '--m','Linewidth',1);
        legend_entries{end+1} = ['V_{plasma,est.Upper} = ',num2str(max(V_plasma_estimate_bounds)),' V'];
        
        plot(min(V_plasma_estimate_bounds).*[1 1], ylim, '--m','Linewidth',1);
        legend_entries{end+1} = ['V_{plasma,est.Lower} = ',num2str(min(V_plasma_estimate_bounds)),' V'];
        
        plot(xlim, I_sat_e_estimate.*[1 1], '-r','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,e,est.}= ',num2str(I_sat_e_estimate),' A'];
        
        plot(xlim, max(I_sat_e_estimate_bounds).*[1 1], '--r','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,e,est.Upper}= ',num2str(max(I_sat_e_estimate_bounds)),' A'];
        
        plot(xlim, min(I_sat_e_estimate_bounds).*[1 1], '--r','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,e,est.Lower}= ',num2str(min(I_sat_e_estimate_bounds)),' A'];
        
        plot(xlim, I_sat_i_estimate.*[1 1], '-b','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,i,est.}= ',num2str(1e3*I_sat_i_estimate),' mA'];
        
        plot(xlim, max(I_sat_i_estimate_bounds).*[1 1], '--b','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,i,est.Upper}= ',num2str(1e3*max(I_sat_i_estimate_bounds)),' mA'];
        
        plot(xlim, min(I_sat_i_estimate_bounds).*[1 1], '--b','Linewidth',1);
        legend_entries{end+1} = ['I_{sat,i,est.Lower}= ',num2str(1e3*min(I_sat_i_estimate_bounds)),' mA'];
        
        ylim([min(I_of_V) max(I_of_V)]);
        [lgd,icons2,plots,txt] = legend(legend_entries,'Location','northwest');
        lgd.FontSize = 14;
        


        T_e_string = ['T_{e,est.} = ',num2str(T_e_estimate,2),' eV ','(',num2str(min(T_e_estimate_bounds),2),' eV',',',num2str(max(T_e_estimate_bounds),2),' eV',')'];
        n_e_string = ['n_{e,est.} = ',num2str(n_e_estimate,'%10.1e'),' m^{-3} ','(',num2str(min(n_e_estimate_bounds),'%10.1e'),' m^{-3}',',',num2str(max(n_e_estimate_bounds),'%10.1e'),' m^{-3}',')'];
        V_plasma_string = ['V_{plasma,est.} = ',num2str(V_plasma_estimate,3),' V ','(',num2str(min(V_plasma_estimate_bounds),3),' V',',',num2str(max(V_plasma_estimate_bounds),3),' V',')'];
        V_float_string = ['V_{float,est.} = ',num2str(V_float_estimate,3),' V ','(',num2str(min(V_float_estimate_bounds),3),' V',',',num2str(max(V_float_estimate_bounds),3),' V',')'];
        text(V_axis(floor(length(V_axis)/6)),I_sat_e_estimate/2,{T_e_string;n_e_string;V_plasma_string;V_float_string},'FontSize',14);
        
        hold off
        
%         V_bars_signal = V_bars(1).*ones(length(pp_stationary.ivt_source.ivs_source_plasma.V_signal),1);
%         I_bars_signal = I_bars(1).*ones(length(pp_stationary.ivt_source.ivs_source_plasma.V_signal),1);
%         
%         figure
%         errorbar(time_axis.*1e6,pp_stationary.ivt_source.ivs_source_plasma.V_signal,V_bars_signal,'k.')
%         title('Bias Voltage Signal');
%         xlabel('time [µs]');
%         ylabel('V_{bias} [V]')
%         figure
%         errorbar(time_axis.*1e6,pp_stationary.ivt_source.ivs_source_plasma.I_signal,I_bars_signal,'k.')
%         title('Probe current Signal');
%         xlabel('time [µs]');
%         ylabel('I_{probe} [A]')
        
%         figure
%         title('Plasma parameters fitted with a linear + exponential model');
%         hold on
%         
%         errorbar(V_axis,I_of_V,(I_ADC_step/2).*ones(1,length(V_axis)),'ko');
%         legend_entries = {'Data points'};
%         
%         plot(V_axis,I_i,'b--');
%         legend_entries{end+1} = 'Ion current';
%         
%         plot(xlim, [0 0], '-k','Linewidth',1);
%         legend_entries{end+1} = 'Zero current line';
%         
%         plot(V_float_fit.*[1 1], ylim, '-g','Linewidth',1);
%         legend_entries{end+1} = ['V_{float,fit} = ',num2str(V_float_fit),' V'];
%         
%         plot(V_axis(knee_ind).*[1 1], ylim, '--k','Linewidth',1);
%         legend_entries{end+1} = 'Inflection threshold';
%         
%         plot(V_plasma_fit.*[1 1], ylim, '-m','Linewidth',1);
%         legend_entries{end+1} = ['V_{plasma,fit} = ',num2str(V_plasma_fit),' V'];
%         
%         plot(xlim, I_sat_e_fit.*[1 1], '-r','Linewidth',1);
%         legend_entries{end+1} = ['I_{sat,e,fit}= ',num2str(I_sat_e_fit),' A'];
%         
%         plot(xlim, I_sat_i_fit.*[1 1], '-b','Linewidth',1);
%         legend_entries{end+1} = ['I_{sat,i,fit}= ',num2str(I_sat_i_fit),' A'];
%         
%         plot(V_axis, I_fit(V_axis), '-k','Linewidth',1);
%         legend_entries{end+1} = 'I_{fit}(V) = I_i + I_{sat,e}*exp[(e/k_BT_e)(V-V_{plasma})]';
%         
%         ylim([min(I_of_V) max(I_of_V)]);
%         [lgd2,icons2,plots,txt] = legend(legend_entries,'Location','northwest');
%         lgd2.FontSize = 14;
%         
%         n_e_fit = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,I_sat_i_fit,T_e_fit);
%         
%         T_e_string = ['T_{e,fit} = ',num2str(T_e_fit,2),' eV'];
%         n_e_string = ['n_{e,fit} = ',num2str(n_e_fit,'%10.1e'),' m^{-3}'];
%         V_plasma_string = ['V_{plasma,fit} = ',num2str(V_plasma_fit,2),' V'];
%         text(min(V_axis)+10,I_sat_e_fit/2,{T_e_string;n_e_string;V_plasma_string},'FontSize',14);
%         
%         
%         hold off
        
        
    case 2
        
%         vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_03_07_Langmuir_Vacuum\f_01KHz\';
%         vacuum_filename = 'meas_0000.h5';
%         vacuum_full_file_path = strcat(vacuum_source_folder,vacuum_filename);

        iv_signals_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_22_Plasma_Profile_100kHz_10V_5_rep_CompleteProfile_X-Drive_On_1500V_Rec-Drive_On_1000V\';
        plasma_filename = 'meas_1040.h5';

        iv_signals_full_file_path = strcat(iv_signals_source_folder,plasma_filename);

        pp_timeseries = PlasmaParamTimeseries(iv_signals_full_file_path);
        nc = 36;
        
        version_selector = 2;
        options.diag_msg_flag = true;
        options.T_e_expected = 10;
        options.T_e_max_allowed = 50;
        options.n_e_expected = 1e18;
        pp_timeseries.extractPlasmaParameters(S_probe,A_ion,version_selector,options,nc,rf);
        
        if ~isempty(pp_timeseries.fitresult_ca{nc,rf})
            
            switch rf
                case 1
                    V_axis = pp_timeseries.ivt_source.rise(nc).V_axis;  
                    I_of_V = pp_timeseries.ivt_source.rise(nc).I_of_V; 
                case 2
                    V_axis = pp_timeseries.ivt_source.fall(nc).V_axis;  
                    I_of_V = pp_timeseries.ivt_source.fall(nc).I_of_V; 
            end
            
            I_of_V_sgolay = pp_timeseries.diagnostics_ca{nc,rf}.I_of_V_sgolay;
            I_of_V_spline = pp_timeseries.diagnostics_ca{nc,rf}.I_of_V_spline;

            pp_I = pp_timeseries.diagnostics_ca{nc,rf}.pp_I;

            p_I = pp_timeseries.diagnostics_ca{nc,rf}.p_I;
            p_DI = pp_timeseries.diagnostics_ca{nc,rf}.p_DI;
            p_D2I = pp_timeseries.diagnostics_ca{nc,rf}.p_D2I;

            logI = pp_timeseries.diagnostics_ca{nc,rf}.logI;
            p_logI = pp_timeseries.diagnostics_ca{nc,rf}.p_logI;
            p_DlogI = pp_timeseries.diagnostics_ca{nc,rf}.p_DlogI;

            I_sat_i_expected = pp_timeseries.diagnostics_ca{nc,rf}.I_sat_i_expected;
            I_sat_e_expected = pp_timeseries.diagnostics_ca{nc,rf}.I_sat_e_expected;
            T_e_expected = pp_timeseries.diagnostics_ca{nc,rf}.T_e_expected;

            V_tf = pp_timeseries.diagnostics_ca{nc,rf}.V_tf;
            I_tf = pp_timeseries.diagnostics_ca{nc,rf}.I_tf;

            I_of_V_fitted = pp_timeseries.I_of_V_fitted_ca{nc,rf};
            fitresult = pp_timeseries.fitresult_ca{nc,rf};

            V_float = pp_timeseries.V_float_a(nc,rf);
            V_plasma = pp_timeseries.V_plasma_a(nc,rf);
            I_sat_i = pp_timeseries.I_sat_i_a(nc,rf);
            I_sat_e = pp_timeseries.I_sat_e_a(nc,rf);
            T_e = pp_timeseries.T_e_a(nc,rf);
            n_e_estimate = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,I_sat_i,T_e)
            
            

            set(0, 'DefaultFigureRenderer', 'opengl');
            
            figure
            movegui(gcf,'northeast');
            axes
            hold on
            
            delta_V = ones(1,length(V_axis)).*pp_timeseries.ivt_source.ivs_source_plasma.V_ADC_step;
            delta_I = ones(1,length(V_axis)).*pp_timeseries.ivt_source.ivs_source_plasma.I_ADC_step;
            errorbar(V_axis,I_of_V,delta_I./2,delta_I./2,delta_V./2,delta_V./2,'o');
            legend_entries = {'Data points'};
            
%             plot(V_axis,I_of_V_sgolay);
%             legend_entries{end+1} = 'I(V)_{sgolay}';

%             plot(V_axis,I_of_V_spline);
%             legend_entries{end+1} = 'I(V)_{spline}';

            %plot(V_axis,I_of_V_fitted,'-k','Linewidth',2);
            %legend_entries{end+1} = 'I(V) = I_{sat,i} + I_{sat,e}*exp((V-V_{plasma})/T_e)';
            hpf = plot(fitresult,V_tf,I_tf);
            delete(hpf(1));
            set(hpf(2),'Color','k');
            set(hpf(2),'Linewidth',2);
            legend_entries{end+1} = 'I(V) = I_{sat,i} + I_{sat,e}*exp((V-V_{plasma})/T_e)';
%             disp(fitresult);
%             I_sat_e_recalculated = -(1/4)*PlasmaPhysics.e*n_e*(1.6)*sqrt(abs(PlasmaPhysics.e)*T_e/PlasmaPhysics.m_e)*S_probe
            
            plot(V_float.*[1 1], ylim, '-g','Linewidth',1);
            legend_entries{end+1} = ['V_{float} = ',num2str(V_float),' V'];
            
            plot(V_plasma.*[1 1], ylim, '-m','Linewidth',1);
            legend_entries{end+1} = ['V_{plasma} = ',num2str(V_plasma),' V'];
            
            plot(xlim, I_sat_i.*[1 1], '-b','Linewidth',1);
            legend_entries{end+1} = ['I_{sat,i}= ',num2str(I_sat_i),' A'];
            
            plot(xlim, I_sat_e.*[1 1], '-r','Linewidth',1);
            legend_entries{end+1} = ['I_{sat,e}= ',num2str(I_sat_e),' A'];
            
            ylim([min(I_of_V) max(I_of_V)]);
            [lgd,icons2,plots,txt] = legend(legend_entries,'Location','northwest');

            hold off
            
%             figure
%             movegui(gcf,'north');
%             axes
%             hold on
%             plot(V_axis,polyval(p_DI,V_axis));
%             hold off
% 
%             figure
%             movegui(gcf,'northeast');
%             axes
%             hold on
%             plot(V_axis,polyval(p_D2I,V_axis));
%             hold off
%             
%             figure
%             movegui(gcf,'east');
%             axes
%             hold on
%             plot(V_axis,log(I_of_V_sgolay - I_sat_i));
%             hold off
            
        end
end



