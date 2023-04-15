

% addpath('..\..\..\General_Purpose');
% addpath('..\..\..\General_Purpose\periodical_cubic_interpolation');
% addpath('..\..\..\Plasma_Physics');
% addpath('..\..\..\Langmuir Analysis V4.0\Classes');



% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_15_ECRH_Langmuir_Profile_10kHz_5V\oscilloscope\';
% plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_22_Plasma_Profile_100kHz_10V_5_rep_CompleteProfile_X-Drive_On_1500V_Rec-Drive_On_1000V\';
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_07_24_Langmuir_ECRH_IV_characteristic_Profile\';


% clearvars -except f_T_e ax_T_e h_T_e_cbar f_n_e ax_n_e h_n_e_cbar f_V_plasma ax_V_plasma h_V_plasma_cbar

load(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
Nc = plasma_ST_profiles.Nt;
load(strcat(plasma_source_folder,'plasma_parameters_multidim_arrays\','T_e_MDa.mat'),'T_e_MDa');
load(strcat(plasma_source_folder,'plasma_parameters_multidim_arrays\','I_sat_i_MDa.mat'),'I_sat_i_MDa');
load(strcat(plasma_source_folder,'plasma_parameters_multidim_arrays\','V_plasma_MDa.mat'),'V_plasma_MDa');

Nx = plasma_ST_profiles.Nx;
Ny = plasma_ST_profiles.Ny;
x_axis = plasma_ST_profiles.x_axis;
y_axis = plasma_ST_profiles.y_axis;

x_pcolor = [x_axis;x_axis(end)+10]-5;
y_pcolor = [y_axis;y_axis(end)+10]-5;

draw_T_e = true;
draw_n_e = true;
draw_V_plasma = true;

keySet = {'stationary','timeseries'};
valueSet = [1 2];
M = containers.Map(keySet,valueSet);

select = 'timeseries'
%select = 'stationary'


switch M(select)
    case 1
        %Stationary mode
        %Average along time_axis (rise)
        T_e_rise = mode(T_e_MDa(:,:,1,:,1),4);
        I_sat_i_rise = mode(I_sat_i_MDa(:,:,1,:,1),4);
        V_plasma_rise = mode(V_plasma_MDa(:,:,1,:,1),4);
        
        %Average along time_axis (fall)
        T_e_fall = mode(T_e_MDa(:,:,1,:,2),4);
        I_sat_i_fall = mode(I_sat_i_MDa(:,:,1,:,2),4);
        V_plasma_fall = mode(V_plasma_MDa(:,:,1,:,2),4);

    case 2
        
        nc_select = 24;
        
        %Timeseries mode
        
        %Toggle movie mode
        movie_mode = false;
        
        %Single frames
        %Select the cycle in the timeseries case
        %nc_select = 9;
        %Pick the selected cycle,average across repetitions (rise)
        T_e_rise      = nanmean(T_e_MDa(:,:,:,nc_select,1),3);
        I_sat_i_rise  = nanmean(I_sat_i_MDa(:,:,:,nc_select,1),3);
        V_plasma_rise = nanmean(V_plasma_MDa(:,:,:,nc_select,1),3);
        %Pick the selected cycle,average across repetitions (fall)
        T_e_fall      = nanmean(T_e_MDa(:,:,:,nc_select,2),3);
        I_sat_i_fall  = nanmean(I_sat_i_MDa(:,:,:,nc_select,2),3);
        V_plasma_fall = nanmean(V_plasma_MDa(:,:,:,nc_select,2),3);
        
        if movie_mode
            %Movies
            %average across repetitions (rise)
            T_e_rise_movie      = nanmean(T_e_MDa(:,:,:,:,1),3);
            I_sat_i_rise_movie  = nanmean(I_sat_i_MDa(:,:,:,:,1),3);
            V_plasma_rise_movie = nanmean(V_plasma_MDa(:,:,:,:,1),3);
            %average across repetitions (fall)
            T_e_fall_movie      = nanmean(T_e_MDa(:,:,:,:,2),3);
            I_sat_i_fall_movie  = nanmean(I_sat_i_MDa(:,:,:,:,2),3);
            V_plasma_fall_movie = nanmean(V_plasma_MDa(:,:,:,:,2),3);
            %Average rise/fall
            clear T_e_movie I_sat_i_movie V_plasma_movie
            for nc = 1:Nc
                T_e_movie{nc} = (T_e_rise_movie(:,:,nc) + T_e_fall_movie(:,:,nc))/2;
                I_sat_i_movie{nc} = (I_sat_i_rise_movie(:,:,nc) + I_sat_i_fall_movie(:,:,nc))/2;
                V_plasma_movie{nc} = (V_plasma_rise_movie(:,:,nc) + V_plasma_fall_movie(:,:,nc))/2;
            end
        end

end
%Average rise/fall
T_e = (T_e_rise + T_e_fall)/2;
I_sat_i = (I_sat_i_rise + I_sat_i_fall)/2;
V_plasma = (V_plasma_rise + V_plasma_fall)/2;


%Calculate electron density n_e[m^-3] from T_e and I_sat_i
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
S_probe = L_tip*2*pi*r_tip;
% Ion mass number
A_ion = 39.948;%Argon
n_e = nan(Nx,Ny);
for nx = 1:Nx
    for ny = 1:Ny
        n_e(nx,ny) = PlasmaPhysics.electron_density_langmuir(S_probe,A_ion,I_sat_i(nx,ny),T_e(nx,ny));
    end
end


%Plasma detector
n_e_min_threshold = 3e18;
for nx = 1:Nx
    for ny = 1:Ny
        if n_e(nx,ny) < n_e_min_threshold || isnan(n_e(nx,ny))
            n_e(nx,ny) = 0;
            T_e(nx,ny) = 0;
            V_plasma(nx,ny) = nan;
        end
    end
end



    %Temperature profile
    if draw_T_e
        n_isolines_T_e = 10;
        if(exist('f_T_e'))
            cla(ax_T_e);
            clf(f_T_e);
        else
            f_T_e = figure;
        end
        figure(f_T_e);
        set(gcf,'Renderer','painters');
        ax_T_e = axes;
        %[~,h_T_e_contourf] = contourf(ax_T_e,x_axis,y_axis,T_e',n_isolines_T_e);
        %set(h_T_e_pcolor,'Linecolor','none');
        h_T_e_imagesc = imagesc(ax_T_e,x_axis,y_axis,T_e');
        set(h_T_e_imagesc,'AlphaData',~isnan(T_e'));
        set(ax_T_e,'YDir','normal');
        caxis(ax_T_e);%, [T_min T_max]);
        set(ax_T_e,'fontsize',20);
        title(ax_T_e,'Electron temperature');
        xlabel(ax_T_e,'X position [mm]');
        ylabel(ax_T_e,'Y position [mm]');
        axis(ax_T_e,'equal');
        colormap(ax_T_e,'hot');
        h_T_e_cbar = colorbar(ax_T_e);
        ylabel(h_T_e_cbar, 'T_e [eV]');
    end

    %Density profile
    if draw_n_e
        %Isolines number
        n_isolines_n_e = 10;

        if(exist('f_n_e'))
            cla(ax_n_e);
            clf(f_n_e);
        else
            f_n_e = figure;
        end
        figure(f_n_e);
        set(gcf,'Renderer','painters');
        ax_n_e = axes;
%         [~,h_n_e_contourf] = contourf(ax_n_e,x_axis,y_axis,n_e',n_isolines_n_e);
%         set(h_n_e_contourf,'Linecolor','none');
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


    %Plasma potential profile
    if draw_V_plasma
        %Isolines number
        n_isolines_V_plasma = 20;
        %Sanity check
        V_min = -20;
        V_max = 20;

        if(exist('f_V_plasma'))
            cla(ax_V_plasma);
            clf(f_V_plasma);
        else
            f_V_plasma = figure;
        end
        figure(f_V_plasma);
        set(gcf,'Renderer','painters');
        ax_V_plasma = axes;
%         [~,h_V_plasma_contourf] = contourf(ax_V_plasma,x_axis,y_axis,V_plasma',n_isolines_V_plasma);
%         set(h_V_plasma_contourf,'Linecolor','none');
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