

if exist('plasma_ST_profiles')
    clearvars -except plasma_ST_profiles
    iv_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';
    v_float_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_Floating_Potential_Profile\';
    st_profile_type = 'stationary';
    repetitions_per_point = 1;
else
    clear
    iv_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';
    v_float_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_Floating_Potential_Profile\';
    st_profile_type = 'stationary';
    repetitions_per_point = 1;
    load(strcat(iv_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles_v7p3.mat'),'plasma_ST_profiles');
    plasma_ST_profiles.buildMDArrays_query;
end

close all

x_axis = plasma_ST_profiles.x_axis;
y_axis = plasma_ST_profiles.y_axis;



% Post-processing and sanity check

%e-Temperature
T_e_rise = plasma_ST_profiles.MD_arrays_ca{1, 2}{1, 1}(:,:,1);
T_e_fall = plasma_ST_profiles.MD_arrays_ca{1, 2}{1, 1}(:,:,2);
T_e = (T_e_rise + T_e_fall) ./2;
clear T_e_rise T_e_fall
T_e(T_e < 0) = nan;
T_e(T_e > 5) = nan;

%e-Density
n_e_rise = plasma_ST_profiles.MD_arrays_ca{2, 2}{1, 1}(:,:,1);
n_e_fall = plasma_ST_profiles.MD_arrays_ca{2, 2}{1, 1}(:,:,2);
n_e = (n_e_rise + n_e_fall) ./2;
clear n_e_rise n_e_fall
n_e(imag(n_e)~=0) = nan;
n_e(n_e < 0) = nan;
n_e(n_e > 1e18) = nan;

%V_plasma
V_plasma_rise = plasma_ST_profiles.MD_arrays_ca{3, 2}{1, 1}(:,:,1);
V_plasma_fall = plasma_ST_profiles.MD_arrays_ca{3, 2}{1, 1}(:,:,2);
V_plasma = (V_plasma_rise + V_plasma_fall) ./2;
clear V_plasma_rise V_plasma_fall
V_plasma(imag(V_plasma)~=0) = nan;
V_plasma(V_plasma < 0) = nan;
V_plasma(V_plasma > 30) = nan;

%V_float
V_float_rise = plasma_ST_profiles.MD_arrays_ca{4, 2}{1, 1}(:,:,1);
V_float_fall = plasma_ST_profiles.MD_arrays_ca{4, 2}{1, 1}(:,:,2);
V_float = (V_float_rise + V_float_fall) ./2;
clear V_float_rise V_float_fall
% V_float(imag(V_plasma)~=0) = nan;
% V_float(V_plasma < 0) = nan;
% V_float(V_plasma > 30) = nan;

%V_plasma
% V_plasma = V_float + 5.2.*T_e;

%Visualization
set(0, 'DefaultFigureRenderer', 'opengl');
draw_T_e = true;
draw_n_e = true;
draw_V_plasma = true;
draw_V_float = true;

%e-Temperature profile
if draw_T_e
    f_T_e = figure;
    ax_T_e = axes;
    h_T_e_imagesc = imagesc(ax_T_e,x_axis,y_axis,T_e');
    set(h_T_e_imagesc,'AlphaData',~isnan(T_e'));
    set(ax_T_e,'YDir','normal');
    caxis(ax_T_e);%, [T_min T_max]);
    set(ax_T_e,'fontsize',20);
    title(ax_T_e,'Electron temperature');
    xlabel(ax_T_e,'X position [mm]');
    ylabel(ax_T_e,'Y position [mm]');
    axis(ax_T_e,'equal');
    colormap(ax_T_e,'autumn');
    h_T_e_cbar = colorbar(ax_T_e);
    ylabel(h_T_e_cbar, 'T_e [eV]');
end

%e-Density profile
if draw_n_e
    %Isolines number
    n_isolines_n_e = 10;
    f_n_e = figure;
    ax_n_e = axes;
    h_n_e_imagesc = imagesc(ax_n_e,x_axis,y_axis,n_e');
    set(h_n_e_imagesc,'AlphaData',~isnan(n_e'));
    set(ax_n_e,'YDir','normal');
    caxis(ax_n_e);%, [n_min n_max]);
    set(ax_n_e,'fontsize',20);
    title(ax_n_e,'Electron density');
    xlabel(ax_n_e,'X position [mm]');
    ylabel(ax_n_e,'Y position [mm]');
    axis(ax_n_e,'equal');
    colormap(ax_n_e,'parula');
    h_n_e_cbar = colorbar(ax_n_e);
    ylabel(h_n_e_cbar, 'n_e [m^{-3}]');
end


%V_plasma profile
if draw_V_plasma
    %Isolines number
    n_isolines_V_plasma = 20;
    f_V_plasma = figure;
    ax_V_plasma = axes;
    h_V_plasma_imagesc = imagesc(ax_V_plasma,x_axis,y_axis,V_plasma');
    set(h_V_plasma_imagesc,'AlphaData',~isnan(V_plasma'));
    set(ax_V_plasma,'YDir','normal');
    caxis(ax_V_plasma);%, [V_min V_max]);
    set(ax_V_plasma,'fontsize',20);
    title(ax_V_plasma,'Plasma Potential');
    xlabel(ax_V_plasma,'X position [mm]');
    ylabel(ax_V_plasma,'Y position [mm]');
    axis(ax_V_plasma,'equal');
    colormap(ax_V_plasma,'jet');
    h_V_plasma_cbar = colorbar(ax_V_plasma);
    ylabel(h_V_plasma_cbar, 'V_{plasma} [V]');
end

%V_float profile
if draw_V_float
    %Isolines number
    n_isolines_V_float = 20;
    f_V_float = figure;
    ax_V_float = axes;
    h_V_float_imagesc = imagesc(ax_V_float,x_axis,y_axis,V_float');
    set(h_V_float_imagesc,'AlphaData',~isnan(V_float'));
    set(ax_V_float,'YDir','normal');
    caxis(ax_V_float);%, [V_min V_max]);
    set(ax_V_float,'fontsize',20);
    title(ax_V_float,'Floating Potential');
    xlabel(ax_V_float,'X position [mm]');
    ylabel(ax_V_float,'Y position [mm]');
    axis(ax_V_float,'equal');
    colormap(ax_V_float,'jet');
    h_V_float_cbar = colorbar(ax_V_float);
    ylabel(h_V_float_cbar, 'V_{float} [V]');
end

