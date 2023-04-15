
addpath('../../../General_Purpose/');
addpath('../../../Plasma_Physics/');
addpath('../../Classes/');

vacuum_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_03_07_Langmuir_Vacuum\f_01KHz\';
vacuum_filename = 'meas_0000.h5';
vacuum_full_file_path = strcat(vacuum_source_folder,vacuum_filename);
% %WARNING: THE DEMPLIFICATION FACTOR DOESNT LOOK RIGHT, I MUST RE-TAKE THE
% %VACUUM MEASUREMENTS

plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_03_07_Langmuir_ECRH_Central_position_triggered\f_01KHz\';
plasma_filename = 'meas_0000.h5';
plasma_full_file_path = strcat(plasma_source_folder,plasma_filename);

pp_stationary = PlasmaParamStationary(plasma_full_file_path,vacuum_full_file_path);
pp_stationary.extractPlasmaParamStationary;




close all

time_axis = pp_stationary.ivt_source.ivs_source_plasma.time_axis;
V_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.V_signal;
I_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.I_signal;

V_signal_vacuum = pp_stationary.ivt_source.ivs_source_vacuum.V_signal;
I_signal_vacuum = pp_stationary.ivt_source.ivs_source_vacuum.I_signal;

V_signal_plasma = pp_stationary.ivt_source.ivs_source_plasma.V_signal;
I_signal_plasma_CM_corrected = pp_stationary.ivt_source.ivs_source_plasma.I_CM_corrected;
I_smooth_exp_fit = pp_stationary.ivt_source.ivs_source_plasma.I_smooth_exp_fit;



figure
set(gcf,'Renderer','painters');
movegui(gcf,'northwest');
subplot(2,1,1);
plot(time_axis,V_signal_plasma);
subplot(2,1,2);
plot(time_axis,I_signal_plasma);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'north');
subplot(2,1,1);
plot(time_axis,V_signal_vacuum);
subplot(2,1,2);
plot(time_axis,I_signal_vacuum);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'west');
subplot(2,1,1);
plot(time_axis,V_signal_plasma);
subplot(2,1,2);
plot(time_axis,I_signal_plasma_CM_corrected);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'center');
subplot(2,1,1);
plot(time_axis,V_signal_plasma);
subplot(2,1,2);
plot(time_axis,I_signal_plasma_CM_corrected - I_signal_vacuum);

figure
set(gcf,'Renderer','painters');
movegui(gcf,'east');
plot(time_axis,I_signal_plasma - I_signal_vacuum);
hold on
plot(time_axis,smooth(I_signal_plasma - I_signal_vacuum,1000,'lowess'));
plot(time_axis,I_smooth_exp_fit);



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



figure
set(gcf,'Renderer','painters');
movegui(gcf,'northeast');
plot(V_axis,I_of_V);
hold
plot(V_axis,I_of_V_sgolay);
plot(V_axis,I_of_V_fitted);
ylim([min(I_of_V) max(I_of_V)]);


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


