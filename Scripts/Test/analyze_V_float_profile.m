
addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\General_Purpose');
addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Plasma_Physics');
addpath('C:\Users\tzf\Desktop\VINETA II Lab\MATLAB\Langmuir Analysis V4.0\Classes');

plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_19_Langmuir_ECRH_IV_characteristic_Profile\';

v_float_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_19_Langmuir_ECRH_Floating_Potential_Profile\';

plasma_st_profiles = PlasmaSTProfiles(plasma_source_folder);
plasma_st_profiles.readMapFile(1);
