
% addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\General_Purpose');
% addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Plasma_Physics');
% addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Langmuir Analysis V4.0\Classes');

%Langmuir probe area and gas used
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_probe = L_tip*2*pi*r_tip;
% Ion mass number
A_ion = 39.948;%Argon

iv_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';
v_float_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_Floating_Potential_Profile\';
st_profile_type = 'stationary';
repetitions_per_point = 1;

plasma_ST_profiles = PlasmaSTProfiles(st_profile_type,iv_source_folder,v_float_source_folder);
plasma_ST_profiles.readMapFile_parallel(repetitions_per_point);

% plasma_ST_profiles.readMapFile(5);
% save(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
% clear plasma_ST_profiles
% load(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');

langmuir_iv_fit_version = 3;
options.diag_msg_flag = false;
options.T_e_expected = 4;
options.T_e_max_allowed = 10;
options.n_e_expected = 1e17;
options.include_external_V_float = true;
% plasma_ST_profiles.extractPlasmaSTProfiles(S_probe,A_ion,langmuir_iv_fit_version,options);
plasma_ST_profiles.extractPlasmaSTProfiles_parallel(S_probe,A_ion,langmuir_iv_fit_version,options);
% plasma_ST_profiles.extractPlasmaSTProfiles_dummy(S_probe,A_ion,langmuir_iv_fit_version,options);
save('-v7.3',strcat(iv_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles_v7p3.mat'),'plasma_ST_profiles');


