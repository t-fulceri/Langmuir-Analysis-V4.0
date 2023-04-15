
addpath('../../../General_Purpose/');
addpath('../../../Plasma_Physics/');
addpath('../../Classes/');

% %Plasma signals
% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_15_ECRH_Langmuir_Profile_10kHz_5V\oscilloscope\';
% % plasma_filename = 'meas_0268.h5'; %Central
% plasma_filename = 'meas_0175.h5';
% plasma_full_file_path = strcat(plasma_source_folder,plasma_filename);
% %Vacuum signals
vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_27_Vacuum_Effects\f_010khz\';
vacuum_filename = 'meas_0000.h5';
vacuum_full_file_path = strcat(vacuum_source_folder,vacuum_filename);
% %WARNING: THE DEMPLIFICATION FACTOR DOESNT LOOK RIGHT, I MUST RE-TAKE THE
% %VACUUM MEASUREMENTS
% 
% if ~exist('pp_timeseries','class')
%     pp_timeseries = PlasmaParamTimeseries(plasma_full_file_path);
%     pp_timeseries.extractPlasmaParameters;
% end
% 
% close all;
% time_axis = pp_timeseries.ivt_source.ivs_source_plasma.time_axis;
% V_signal = pp_timeseries.ivt_source.ivs_source_plasma.V_signal;
% V_ADC_step = pp_timeseries.ivt_source.ivs_source_plasma.V_ADC_step;
% I_signal = pp_timeseries.ivt_source.ivs_source_plasma.I_signal;
% I_smooth = pp_timeseries.ivt_source.ivs_source_plasma.I_smooth;
% I_smooth_exp_fit = pp_timeseries.ivt_source.ivs_source_plasma.I_smooth_exp_fit;
% I_CM_corrected = pp_timeseries.ivt_source.ivs_source_plasma.I_CM_corrected;
% time_axis_rise = pp_timeseries.ivt_source.total_time_axis_rise;
% time_axis_fall = pp_timeseries.ivt_source.total_time_axis_fall;


plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_15_ECRH_Langmuir_Profile_10kHz_5V\oscilloscope\';
plasma_filename = 'meas_0020.h5'; %Central
%plasma_filename = 'meas_0175.h5';
plasma_full_file_path = strcat(plasma_source_folder,plasma_filename);

% plasma_ST_profiles = PlasmaSTProfiles(plasma_source_folder);
% plasma_ST_profiles.readMapFile;
% disp(plasma_ST_profiles);
% %plasma_ST_profiles.extractPlasmaSTProfiles;
% save(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
% clear plasma_ST_profiles

% load(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
% plasma_ST_profiles.buildMDArrays;

%Langmuir probe area and gas used
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_probe = L_tip*2*pi*r_tip;
% Ion mass number
A_ion = 39.948;%Argon

pp_stationary = PlasmaParamStationary(plasma_full_file_path,vacuum_full_file_path);
pp_stationary.extractPlasmaParameters(S_probe,A_ion,2);

rf = 2;
switch rf
    case 1
        V_axis = pp_stationary.V_axis_rise;
        I_of_V = pp_stationary.I_of_V_rise;
    case 2
        V_axis = pp_stationary.V_axis_fall;
        I_of_V = pp_stationary.I_of_V_fall;
end
I_of_V_sgolay = pp_stationary.I_of_V_sgolay_ca{rf};
I_of_V_fitted = pp_stationary.I_of_V_fitted_ca{rf};


close all
plot(V_axis,I_of_V);
hold
plot(V_axis,I_of_V_sgolay);
plot(V_axis,I_of_V_fitted);


T_e = pp_stationary.T_e(rf)
I_sat_i = pp_stationary.I_sat_i(rf)

%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S = L_tip*2*pi*r_tip;
% Ion mass number
A = 39.948;%Argon

pp = PlasmaPhysics;
pp.electron_density_langmuir(S,A,I_sat_i,T_e);


