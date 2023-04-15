% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_06_05_Langmuir_ECRH_Power_Cycle_Characterization_New_Meas_Config\f_500hz_plasma\';
% plasma_filename = 'meas_0001.h5';
% plasma_full_file_path = strcat(plasma_source_folder,plasma_filename);
% pp_stationary = PlasmaParamStationary(plasma_full_file_path);

iv_signals_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';
vf_signal_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_Floating_Potential_Profile\';
plasma_filename = 'meas_0200.h5';

iv_signals_full_file_path = strcat(iv_signals_source_folder,plasma_filename);
vf_signal_full_file_path = strcat(vf_signal_source_folder,plasma_filename);
pp_stationary = PlasmaParamStationary(iv_signals_full_file_path,vf_signal_full_file_path);

close all
figure
hold on
plot(pp_stationary.V_axis_rise,pp_stationary.I_of_V_rise,'or');
plot(pp_stationary.V_axis_fall,pp_stationary.I_of_V_fall,'ok');
hold off