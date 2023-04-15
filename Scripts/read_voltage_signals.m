

%plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_06_04_Post_Isolator_Voltage_Offset\';
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_19_Langmuir_ECRH_Floating_Potential_Profile\';
[file_name,file_folder] =  uigetfile([plasma_source_folder,'*.h5'],plasma_source_folder);
full_file_path = [file_folder,file_name];

v_signal = VSignal(full_file_path);