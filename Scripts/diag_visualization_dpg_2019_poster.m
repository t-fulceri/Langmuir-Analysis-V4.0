close all;

%% Load Langmuir and B-dot data

% %X-Drive 1500V and Rex-Drive 1000V folders
langmuir_data_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_22_Plasma_Profile_100kHz_10V_5_rep_CompleteProfile_X-Drive_On_1500V_Rec-Drive_On_1000V\';
b_dot_data_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_B-dot_CompleteProfile_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\01\';

%ECRH ONLY with dummy B_dot
% langmuir_data_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\Reconstructed_Profiles\2018_08_21_Plasma_Profile_01kHz_10V_1_rep_CurrentMonitor_CompleteProfile_ECRH_Only_DeamplifiedBy200\';
% b_dot_data_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2018_08_24_B-dot_CompleteProfile_E-Guns_On_X-Drive_On_1500V_Rec-Drive_On_1000V\01\';

%% Load Reconnection monitor data
if ~exist('rec_data')
    rec_data = load([b_dot_data_folder,'recData.mat']);
end
%Load the timeseries data structure
if ~exist('iv_timeseries_cell_array')
    load([langmuir_data_folder, 'iv_timeseries_cell_array']);
    iv_timeseries = iv_timeseries_cell_array{1,1};
end

%% Load the profiles

load([langmuir_data_folder, 'profiles_to_plot']);

%% Information on drawing axes

%Draw range
draw_range_x = [-150,+150];
draw_range_y = [-150,+150];

%Langmuir origin WRT axis (from laboratory)
x_L_0 = -148;
y_L_0 = -115;

%B-dot origin WRT axis (from laboratory)
x_B_0 = -138;
y_B_0 = -75;

%Scan range
scan_range_x = 0:10:220;
scan_range_y = 0:10:170;

%Langmuir scan range
scan_range_L_x = scan_range_x + x_L_0;
scan_range_L_y = scan_range_y + y_L_0;

%B-dot scan range
scan_range_B_x = scan_range_x + x_B_0;
scan_range_B_y = scan_range_y + y_B_0;



%% Physical Constants

%Electron_charge [C]
e = -1.6021766208 * 1e-19;
%Boltzmann constant [J/K]
k_B = 1.38064852 * 1e-23;
%Atomic Mass Unit [kg]
AMU = 1.660539040*1e-27;

%% Langmuir probe surface and ion mass number
%Probe tip area
L_tip = 4.00*1e-3; %4.00 mm length
r_tip = 0.25*1e-3; %0.25 mm radius
probe_tip_area = L_tip*2*pi*r_tip;
% Ion mass number
ion_mass_number = 39.948;%Argon

%% Convert ion saturation current into electron density

%Convert electron temperature(s) to SI units (Kelvin)
T_e_XYC_SI = T_e_XYC.* ((-e)/k_B);
%Convert ion mass from AMU to kg
m_i = ion_mass_number*AMU;
%Calculate electron charge density/ies [C*(m^-3)]
rho_e_XYC = (I_sat_i_XYC./probe_tip_area) .* sqrt(2*pi*m_i./(k_B.*T_e_XYC_SI));
%Calculate electron density [m^-3]
n_e_XYC = rho_e_XYC./e;
%----------------------------------------------------------------

%% Pressure profile 
p_XYC = k_B.*n_e_XYC.*T_e_XYC;

%% For Ion Saturation Current: Reverse sign and convert from A to mA
I_sat_i_XYC = -1e3.*I_sat_i_XYC;

%% T_e, V_n_e, and V_plasma scales
T_e_scale = [5 15];
V_plasma_scale = [min(min(min(V_plasma_XYC))),max(max(max(V_plasma_XYC)))];
n_e_scale = [0 1e20];

%% Smooth the data a little bit...
space_smoothing = [3,3,1];
time_smoothing = [1,1,3];
spacetime_smoothing = 3;
%T_e_XYC = smooth3(T_e_XYC,'box',space_smoothing);
% I_sat_i_XYC = smooth3(I_sat_i_XYC,'box',space_smoothing);
V_plasma_XYC = smooth3(V_plasma_XYC,'box',space_smoothing);

%n_e_XYC = smooth3(n_e_XYC,'box',space_smoothing);
p_XYC = smooth3(p_XYC,'box',space_smoothing);

%Axial current profile
j_z_XYC = smooth3(rec_data.j,'box',space_smoothing);
%Axial vector potential profile
A_z_XYC = smooth3(rec_data.AP,'box',space_smoothing);
%Axial electric field profile
E_z_XYC = diff(rec_data.AP,200,3);


%% Set isolines number
n_isolines_T_e = 10;
n_isolines_n_e = 10;
n_isolines_V_plasma = 10;
n_isolines_Az = 30;
n_isolines_jz = 10;
csp_smoothing_param = 0.5;


%% Timing information
%snapshots = [4 ,8, 12, 16];
snapshots = 1:iv_timeseries.NC_tot;

Ns_tot = length(snapshots);
ns = 18;

% %Not exactly the same but ok for testing
% time_axis = 1e6.*(iv_timeseries.total_time_axis_rise);

time_axis = 1e6.*(iv_timeseries.total_time_axis_rise + iv_timeseries.total_time_axis_fall)./2;

%Extract the time elapsed since LABVIEW trigger [seconds]
real_timestamp = time_axis(ns)*1e-6 + 20e-3;

%Look for corresponding index in B-dot data timeseries
%plasData
[~,plasData_time_index] = min(abs(real_timestamp-rec_data.time.plasData));
%vacData
[~,vacData_time_index] = min(abs(real_timestamp-rec_data.time.vacData));
%Reconnection Control Center (rcc)
[~,rcc_time_index] = min(abs(real_timestamp-rec_data.time.rcc));

%Record value of X-Drive at selected time_index
mean_U_timeseries = 2.*mean(rec_data.rccPlas(15).data,2);
x_drive_current_timeseries = (rec_data.rccPlas(15).const)*cumsum(mean_U_timeseries).*1e-10;
x_drive_current = x_drive_current_timeseries(rcc_time_index);


%% Plot selction switches
%Reconnection monitor
cathode_curent_plot = false;
xdrive_and_recdrive_plot = false;
%Langmuir diagnostic
T_e_plot = true;
n_e_plot = true;
V_plasma_plot = true;
p_and_Exy_plot = true;

B_field_lines = true;

T_e_plot_ECRH = false;
n_e_plot_ECRH = false;

%B-dot diagnostic
A_z_plot = false;
A_z_and_j_z_plot = true;
E_z_plot = false;

n_e_and_A_z_plot = false;
T_e_and_j_z_plot = false;

%% Movie flag
movie_flag = false;

%% Clear before drawing
close all
clear xlim ylim

%% Common settings
font_size_poster = 20;
current_traces_t_limits = [0 300];
e_gun_traces_I_limits = [0 1500];
xd_recd_traces_I_limits = [-1000 2000];
current_fig_width = 800;
current_fig_height = 400;

x_axis_label = 'X position [mm]';
y_axis_label = 'Y position [mm]';
profiles_fig_x = 100;
profiles_fig_y = 100;
profiles_fig_width = 800;
profiles_fig_height = 600;
Nx = length(scan_range_L_x);
Ny = length(scan_range_L_y);
draw_limits_L_x = [min(scan_range_L_x) max(scan_range_L_x)];
draw_limits_L_y = [min(scan_range_L_y) max(scan_range_L_y)];



%% cathode_curent_plot
if cathode_curent_plot
    %-------------CATHODE CURRENT--------------
    figure
    set(gcf,'Renderer','painters');
    movegui(gcf,'center');

    cathode_current = mean(rec_data.rccPlas(8).data,2) - mean(rec_data.rccVac(8).data,2);
    
    plot((rec_data.time.rcc - 20e-3).*1e6,-cathode_current,'LineWidth',2);
    hold on
    xlim(current_traces_t_limits);
    ylim(e_gun_traces_I_limits);
    plot([time_axis(ns) time_axis(ns)], ylim, '-r','LineWidth',2);
    hold off
    
    set(gca,'fontsize',font_size_poster);

    title('E-gun Cathode Current Timetrace');
    xlabel('Time [탎]');
    ylabel('Current [A]');
    %------------------------------------------
    set(gcf, 'Position', [100, 100, 100 + current_fig_width, 100 + current_fig_height]);
end

%% xdrive_and_recdrive_plot
if xdrive_and_recdrive_plot
    %-------------X_DRIVE AND REC_DRIVE CURRENT--------------
    figure
    set(gcf,'Renderer','painters');
    movegui(gcf,'center');

    x_drive_current = -rec_data.rccVac(15).data(:,1);
    rec_drive_current = -rec_data.rccVac(16).data(:,1);
    
    plot((rec_data.time.rcc - 20e-3).*1e6,x_drive_current,'LineWidth',2,'Color',[1 0.5 0]);
    hold on
    plot((rec_data.time.rcc - 20e-3).*1e6,rec_drive_current,'LineWidth',2,'Color','blue');
    xlim(current_traces_t_limits);
    ylim(xd_recd_traces_I_limits);
    plot([time_axis(ns) time_axis(ns)], ylim, '-r','LineWidth',2);
    set(gca,'fontsize',font_size_poster);
    
    
    title('X-drive and Reconnection-drive Current');
    %title(['\color[rgb]{1 .5 0}X-drive', '\color{black} and ', '\color{blue}Reconnection-drive','\color{black} currents'])
    xlabel('Time [탎]');
    ylabel('Current [A]');
    hold off
    %------------------------------------------
    legend('X-drive current','Rec-drive current');
    set(gcf, 'Position', [100, 100, 100 + current_fig_width, 100 + current_fig_height]);
end


%% T_e_plot
if T_e_plot
    T_e_f = figure;
    T_e_ax = axes;
    %Avoid graphical errors
    set(T_e_f,'Renderer','painters');
    [~,T_e_cf] = contourf(scan_range_L_x,scan_range_L_y,(T_e_XYC(:,:,snapshots(ns)))',n_isolines_T_e);
    set(T_e_cf,'LineColor','none');
    set(T_e_ax,'fontsize',font_size_poster);
    %title(strcat('Electron temperature @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title(T_e_ax,'Electron temperature');
    xlabel(T_e_ax,x_axis_label);
    ylabel(T_e_ax,y_axis_label);
    colormap(T_e_ax,hot);
    %caxis(T_e_scale);
    T_e_cb = colorbar(T_e_ax);
    ylabel(T_e_cb, 'T_e [eV]');
    hold on;
    plot(T_e_ax,[0 0], ylim, '-k','Linewidth',1);
    plot(T_e_ax,xlim, [0 0], '-k','Linewidth',1);
    
    if B_field_lines == true
        [~,A_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
        set(A_z_cb,'LineColor','k','LineWidth',1.5);
        xlim(draw_limits_L_x);
        ylim(draw_limits_L_y);
        title(T_e_ax,'Electron temperature and in-plane B-field lines');
    end
    
    set(T_e_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    pbaspect(T_e_ax,[Nx Ny 1]);
end

%% n_e_plot
if n_e_plot
    
    %Conversion (Olaf Aesthetics)
    n_e_XYC = n_e_XYC./1e19;
    
    n_e_f = figure;
    n_e_ax = axes;
    %Avoid graphical errors
    set(n_e_f,'Renderer','painters');
    [~,n_e_cf] = contourf(scan_range_L_x,scan_range_L_y,(n_e_XYC(:,:,snapshots(ns)))',n_isolines_n_e);
    set(n_e_cf,'LineColor','none');
    set(n_e_ax,'fontsize',font_size_poster);
    %title(strcat('Electron density @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Electron density');
    xlabel(n_e_ax,x_axis_label);
    ylabel(n_e_ax,y_axis_label);
    colormap(n_e_ax,parula);
    %caxis(n_e_scale);
    n_e_cb = colorbar(n_e_ax);
    ylabel(n_e_cb, 'n_e [10^{19} m^{^-3}]');
    hold on;
    plot([0 0], ylim, '-k');
    plot(xlim, [0 0], '-k');
    
    if B_field_lines == true
        [~,A_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
        set(A_z_cb,'LineColor','k','LineWidth',1.5);
        xlim(draw_limits_L_x);
        ylim(draw_limits_L_y);
        title(n_e_ax,'Electron density and in-plane B-field lines');
    end
    
    set(n_e_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    pbaspect(n_e_ax,[Nx Ny 1]);
end

%% V_plasma_plot
if V_plasma_plot
    V_plasma_f = figure;
    V_plasma_ax = axes;
    %Avoid graphical errors
    set(V_plasma_f,'Renderer','painters');
    [~,V_plasma_cf] = contourf(scan_range_L_x,scan_range_L_y,(V_plasma_XYC(:,:,snapshots(ns)))',n_isolines_V_plasma);
    set(V_plasma_cf,'LineColor','k');
    set(V_plasma_ax,'fontsize',font_size_poster);
    %title(strcat('Electron density @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Plasma Potential');
    xlabel(V_plasma_ax,x_axis_label);
    ylabel(V_plasma_ax,y_axis_label);
    cmap1 = colormap(V_plasma_ax,bone);
    %caxis(V_plasma_scale);
    V_plasma_cb = colorbar(V_plasma_ax);
    ylabel(V_plasma_cb, 'V_plasma [V]');
    hold on;
    plot([0 0], ylim, '-k');
    plot(xlim, [0 0], '-k');
    
    if B_field_lines == true
        [~,A_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
        set(A_z_cb,'LineColor','k','LineWidth',1.5);
        xlim(draw_limits_L_x);
        ylim(draw_limits_L_y);
    end
    
    set(V_plasma_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    pbaspect(V_plasma_ax,[Nx Ny 1]);
end

%% p_and_Exy_plot
if p_and_Exy_plot
    p_and_Exy_f = figure;
    p_and_Exy_ax = axes;
    
    %Conversion (Olaf Aesthetics)
    p_XYC = p_XYC./1e-3;
    
    %Avoid graphical errors
    set(p_and_Exy_f,'Renderer','painters');
    [~,p_and_Exy_cf] = contourf(scan_range_L_x,scan_range_L_y,(p_XYC(:,:,snapshots(ns)))');
    hold(p_and_Exy_ax,'on');
    [neg_E_x,neg_E_y] = gradient(V_plasma_XYC(:,:,snapshots(ns))');
    E_x = - neg_E_x;
    E_y = - neg_E_y;
    hq = quiver(scan_range_L_x,scan_range_L_y,E_x,E_y,'color',[1 1 1],'LineWidth',2);
    set(gca,'fontsize',font_size_poster);
    %title(strcat('Pressure and in-plane Electric Field @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Pressure, in-plane E-field, in-plane B-field lines');
    xlabel(p_and_Exy_ax,x_axis_label);
    ylabel(p_and_Exy_ax,y_axis_label);
    colormap(p_and_Exy_ax,jet);
    %caxis(p_scale);
    p_and_Exy_cb = colorbar(gca);
    ylabel(p_and_Exy_cb, 'p [mPa]');
    hold on;
    plot(p_and_Exy_ax,[0 0], ylim, '-k','Linewidth',1);
    plot(p_and_Exy_ax,xlim, [0 0], '-k','Linewidth',1);
    
    if B_field_lines == true
        [~,A_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
        set(A_z_cb,'LineColor','k','LineWidth',1.5);
        xlim(draw_limits_L_x);
        ylim(draw_limits_L_y);
    end
    
    set(p_and_Exy_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    pbaspect(p_and_Exy_ax,[Nx Ny 1]);
end

%% T_e_plot_ECRH
if T_e_plot_ECRH
    figure
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    movegui(gcf,'north');
    %n_e_ECRH_XY = nan(NX,NY);
    T_e_ECRH_XY = nanmean(T_e_XYC,3);
    ind_to_kill = find(T_e_ECRH_XY>15);
    T_e_ECRH_XY(ind_to_kill) = nan;
    [~,T_e_cb] = contourf(scan_range_L_x,scan_range_L_y,(T_e_ECRH_XY)',n_isolines_T_e);
    set(gca,'fontsize',20);
    set(T_e_cb,'LineColor','none');

    %title(strcat('Electron density @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Electron temperature of ECRH Plasma');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,hot);
    %caxis(T_e_scale);
    T_e_cb = colorbar(gca);
    ylabel(T_e_cb, 'T_e [eV]');
    xlim(draw_range_x);
    ylim(draw_range_y);
    set(gcf, 'Position', [100, 100, 900, 900]);

    hold on;

    plot([0 0], ylim, '-k');
    plot(xlim, [0 0], '-k');
end

%% n_e_plot_ECRH
if n_e_plot_ECRH
    figure
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    movegui(gcf,'north');
    %n_e_ECRH_XY = nan(NX,NY);
    n_e_ECRH_XY = nanmean(n_e_XYC,3);
    [~,T_e_cb] = contourf(scan_range_L_x,scan_range_L_y,(n_e_ECRH_XY)',n_isolines_n_e);
    set(gca,'fontsize',20);
    set(T_e_cb,'LineColor','none');

    %title(strcat('Electron density @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Electron density of ECRH plasma');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,parula);
    %caxis(n_e_scale);
    T_e_cb = colorbar(gca);
    ylabel(T_e_cb, 'n_e [m^{^-3}]');
    xlim(draw_range_x);
    ylim(draw_range_y);
    set(gcf, 'Position', [100, 100, 900, 900]);

    hold on;

    plot([0 0], ylim, '-k');
    plot(xlim, [0 0], '-k');
end

%% A_z_plot
if A_z_plot
    A_z_f = figure;
    A_z_ax = axes;
    %Avoid graphical errors
    set(A_z_f,'Renderer','painters');
    
    [C,A_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
    set(A_z_ax,'fontsize',font_size_poster);
    set(A_z_cb,'LineColor','k','LineWidth',1.5);

    % title(strcat('Electron Density and Axial Vector Potential @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Magnetic field lines');
    xlabel(x_axis_label);
    ylabel(y_axis_label);
    set(A_z_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    hold on;
    plot(A_z_ax,[0 0], ylim, '-k','Linewidth',1);
    plot(A_z_ax,xlim, [0 0], '-k','Linewidth',1);
    pbaspect(A_z_ax,[Nx Ny 1]);
end

%% A_z_and_j_z_plot
if A_z_and_j_z_plot
    A_z_and_j_z_f = figure;
    A_z_and_j_z_ax = axes;
    %Avoid graphical errors
    set(A_z_and_j_z_f,'Renderer','painters');
    movegui(A_z_and_j_z_f,'southeast')

    [~,A_z_and_j_z_cb] = contourf(scan_range_B_x,scan_range_B_y,csp(j_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_n_e);
    set(A_z_and_j_z_cb,'LineColor','none');

    colormap(gca,jet);
    A_z_and_j_z_cb = colorbar(gca);
    ylabel(A_z_and_j_z_cb, 'j_z [A/m^2]');

    hold on;

    [C,A_z_and_j_z_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_Az);
    set(gca,'fontsize',font_size_poster);
    set(A_z_and_j_z_cb,'LineColor','k','LineWidth',1.5);

    % title(strcat('Electron Density and Axial Vector Potential @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Axial current density and in-plane B-field lines');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
%     pbaspect([1 1 1]);
%     xlim(draw_range_x);
%     ylim(draw_range_y);
    set(A_z_and_j_z_f, 'Position', [profiles_fig_x, profiles_fig_y, profiles_fig_x + profiles_fig_width, profiles_fig_y + profiles_fig_height]);
    pbaspect(A_z_and_j_z_ax,[Nx Ny 1]);

    hold on;

    plot([0 0], ylim, '-k','Linewidth',1);
    plot(xlim, [0 0], '-k','Linewidth',1);
end

%% E_z_plot
if E_z_plot
    clear h C;
    figure;
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    movegui(gcf,'southeast')
    
    %X-Drive cross-section
    X_drive_radius = 4;

    th = 0:pi/50:2*pi;
    xunit = X_drive_radius * cos(th) +75;
    yunit = X_drive_radius * sin(th) -130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on
    xunit = X_drive_radius * cos(th) -75;
    yunit = X_drive_radius * sin(th) +130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on

    %[~,h] = contourf(scan_range_B_x,scan_range_B_y,csp(E_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_n_e);
    [~,T_e_cb] = contourf(scan_range_B_x,scan_range_B_y,E_z_XYC(:,:,plasData_time_index),n_isolines_n_e);
    set(T_e_cb,'LineColor','none');
    set(gca,'fontsize',font_size_poster);
    colormap(gca,parula);
    T_e_cb = colorbar(gca);
    ylabel(T_e_cb, 'E_z [V/m]');

    title('Axial Electric Field');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    xlim(draw_range_x);
    ylim(draw_range_y);
    set(gcf, 'Position', [1000, 100, 900, 900]);

    hold on;

    plot([0 0], ylim, '-k','Linewidth',1);
    plot(xlim, [0 0], '-k','Linewidth',1);
end

%% n_e_and_A_z_plot
if n_e_and_A_z_plot
    clear h C;
    figure;
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    movegui(gcf,'southeast')
    
    %X-Drive cross-section
    X_drive_radius = 4;

    th = 0:pi/50:2*pi;
    xunit = X_drive_radius * cos(th) +75;
    yunit = X_drive_radius * sin(th) -130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on
    xunit = X_drive_radius * cos(th) -75;
    yunit = X_drive_radius * sin(th) +130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on

    [~,T_e_cb] = contourf(scan_range_L_x,scan_range_L_y,(n_e_XYC(:,:,snapshots(ns)))',n_isolines_n_e);
    set(T_e_cb,'LineColor','none');

    colormap(gca,parula);
    caxis(n_e_scale);
    T_e_cb = colorbar(gca);
    ylabel(T_e_cb, 'n_e [m^{-3}]');

    hold on;

    [C,T_e_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_jz);
    set(gca,'fontsize',20);
    set(T_e_cb,'LineColor','k','LineWidth',1.5);

    % title(strcat('Electron Density and Axial Vector Potential @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Electron Density and Axial Vector Potential');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    xlim(draw_range_x);
    ylim(draw_range_y);
    set(gcf, 'Position', [100, 100, 900, 900]);

    hold on;

    plot([0 0], ylim, '-k','Linewidth',1);
    plot(xlim, [0 0], '-k','Linewidth',1);
end

%% T_e_and_j_z_plot
if T_e_and_j_z_plot
    
    clear h C;
    figure;
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    movegui(gcf,'south')
    
    %X-Drive cross-section
    X_drive_radius = 4;

    th = 0:pi/50:2*pi;
    xunit = X_drive_radius * cos(th) +75;
    yunit = X_drive_radius * sin(th) -130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on
    xunit = X_drive_radius * cos(th) -75;
    yunit = X_drive_radius * sin(th) +130;
    plot(xunit, yunit,'r','Linewidth',1.5);
    hold on

    [~,T_e_cb] = contourf(scan_range_L_x,scan_range_L_y,(T_e_XYC(:,:,snapshots(ns)))',n_isolines_T_e);
    set(T_e_cb,'LineColor','none');

    set(gca,'fontsize',20);
    %title(strcat('Electron temperature @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
    title('Electron temperature');
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    colormap(gca,hot);
    caxis(T_e_scale);
    T_e_cb = colorbar(gca);
    ylabel(T_e_cb, 'T_e [eV]');

    hold on;
    clear h C;
    [C,T_e_cb] = contour(scan_range_B_x,scan_range_B_y,csp(j_z_XYC(:,:,plasData_time_index),csp_smoothing_param),n_isolines_jz);
    set(gca,'fontsize',20);
    set(T_e_cb,'LineColor','black','LineWidth',1.5);

    title(strcat('Electron Temperature and Axial Current Density'));
    xlabel('X position [mm]');
    ylabel('Y position [mm]');
    pbaspect([1 1 1]);
    xlim(draw_range_x);
    ylim(draw_range_y);
    set(gcf, 'Position', [500, 100, 900, 900]);
    
    hold on;

    plot([0 0], ylim, '-k','Linewidth',1);
    plot(xlim, [0 0], '-k','Linewidth',1);
end

%% Movie production
if movie_flag
    f = figure;
    %Avoid graphical errors
    set(gcf,'Renderer','painters');
    %movegui(gcf,'center');
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf, 'Position', [0, 0, 1600, 1200]);

    frame_rep = 1;
    %Record video
    frames = cell(length(time_axis)*frame_rep,1);
    
    ts_real_begin = 20e-3 + 100e-6;
    ts_real_end = 20e-3 + 250e-6;

    %Look for corresponding index in B-dot data timeseries
    [~,ts_bdot_begin] = min(abs(ts_real_begin-rec_data.time.plasData));
    [~,ts_bdot_end] = min(abs(ts_real_end-rec_data.time.plasData));
    
    for ts_bdot = ts_bdot_begin:ts_bdot_end
        
        ts_real = (rec_data.time.plasData(ts_bdot) - 20e-3) * 1e+6;
        [~,ts_langmuir] = min(abs(ts_real-time_axis));

        subplot(3,2,[1,3]);
        %Avoid graphical errors
        set(gcf,'Renderer','painters');
        %movegui(gcf,'southeast')

        %X-Drive cross-section
        X_drive_radius = 4;

        th = 0:pi/50:2*pi;
        xunit = X_drive_radius * cos(th) +75;
        yunit = X_drive_radius * sin(th) -130;
        plot(xunit, yunit,'r','Linewidth',1.5);
        hold on
        xunit = X_drive_radius * cos(th) -75;
        yunit = X_drive_radius * sin(th) +130;
        plot(xunit, yunit,'r','Linewidth',1.5);
        hold on

        [~,T_e_cb] = contourf(scan_range_B_x,scan_range_B_y,csp(j_z_XYC(:,:,ts_bdot),csp_smoothing_param),n_isolines_n_e);
        set(T_e_cb,'LineColor','none');

        colormap(gca,jet);
        T_e_cb = colorbar(gca);
        ylabel(T_e_cb, 'j_z [A/m^2]');

        hold on;

        [C,T_e_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,ts_bdot),csp_smoothing_param),n_isolines_Az);
        set(gca,'fontsize',10);
        set(T_e_cb,'LineColor','k','LineWidth',1.5);
        caxis([-500 2500]);

        % title(strcat('Electron Density and Axial Vector Potential @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
        title('Axial Vector Potential and Axial Current');
        xlabel('X position [mm]');
        ylabel('Y position [mm]');
        pbaspect([1 1 1]);
        xlim(draw_range_x);
        ylim(draw_range_y);
        %set(gcf, 'Position', [100, 100, 900, 900]);

        hold on;

        plot([0 0], ylim, '-k','Linewidth',1);
        plot(xlim, [0 0], '-k','Linewidth',1);
%         %Avoid graphical errors
%         set(gcf,'Renderer','painters');
%         movegui(gcf,'northwest');
% 
%         [~,h] = contourf(x_range - langmuir_offset_x,y_range - langmuir_offset_y,(T_e_XYC(:,:,ts))',n_isolines_T_e);
%         set(h,'LineColor','none');
% 
%         set(gca,'fontsize',20);
%         %title(strcat('Electron temperature @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
%         title('Electron temperature');
%         xlabel('X position [mm]');
%         ylabel('Y position [mm]');
%         pbaspect([1 length(y_range)/length(x_range) 1]);
%         colormap(gca,hot);
%         caxis(T_e_scale);
%         h = colorbar(gca);
%         ylabel(h, 'T_e [eV]');
% 
%         hold on;

        plot([0 0], ylim, '-k');
        plot(xlim, [0 0], '-k');

        subplot(3,2,[2,4]);
%         %Avoid graphical errors
%         set(gcf,'Renderer','painters');
%         movegui(gcf,'north');
%         [~,h] = contourf(x_range - langmuir_offset_x,y_range - langmuir_offset_y,(n_e_XYC(:,:,ts))',n_isolines_T_e);
%         set(gca,'fontsize',20);
%         %set(h,'LineColor','none');
% 
%         %title(strcat('Electron density @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
%         title('Electron density');
%         xlabel('X position [mm]');
%         ylabel('Y position [mm]');
%         pbaspect([1 length(y_range)/length(x_range) 1]);
%         colormap(gca,parula);
%         caxis(n_e_scale);
%         h = colorbar(gca);
%         ylabel(h, 'n_e [m^{^-3}]');
% 
%         hold on;
% 
%         plot([0 0], ylim, '-k');
%         plot(xlim, [0 0], '-k');

        set(gcf,'Renderer','painters');

        %X-Drive cross-section
        X_drive_radius = 4;

        th = 0:pi/50:2*pi;
        xunit = X_drive_radius * cos(th) +75;
        yunit = X_drive_radius * sin(th) -130;
        plot(xunit, yunit,'r','Linewidth',1.5);
        hold on
        xunit = X_drive_radius * cos(th) -75;
        yunit = X_drive_radius * sin(th) +130;
        plot(xunit, yunit,'r','Linewidth',1.5);
        hold on

        [~,T_e_cb] = contourf(scan_range_L_x,scan_range_L_y,(T_e_XYC(:,:,ts_langmuir))',n_isolines_T_e);
        set(T_e_cb,'LineColor','none');

        set(gca,'fontsize',10);
        %title(strcat('Electron temperature @ t = ',num2str(time_axis(snapshots(ns))),' 탎'));
        title('Electron temperature');
        xlabel('X position [mm]');
        ylabel('Y position [mm]');
        pbaspect([1 1 1]);
        colormap(gca,hot);
        caxis(T_e_scale);
        T_e_cb = colorbar(gca);
        ylabel(T_e_cb, 'T_e [eV]');

        hold on;
        clear h C;
        [C,T_e_cb] = contour(scan_range_B_x,scan_range_B_y,csp(A_z_XYC(:,:,ts_bdot),csp_smoothing_param),n_isolines_Az);
        set(gca,'fontsize',10);
        set(T_e_cb,'LineColor','black','LineWidth',1.5);

        title(strcat('Electron Temperature and Axial Vector Potential'));
        xlabel('X position [mm]');
        ylabel('Y position [mm]');
        pbaspect([1 1 1]);
        xlim(draw_range_x);
        ylim(draw_range_y);
        %set(gcf, 'Position', [500, 100, 900, 900]);

        hold on;

        plot([0 0], ylim, '-k','Linewidth',1);
        plot(xlim, [0 0], '-k','Linewidth',1);

        subplot(3,2,5:6);
        set(gcf,'Renderer','painters');

        x_drive_current = -rec_data.rccVac(15).data(:,1);
        rec_drive_current = -rec_data.rccVac(16).data(:,1);

        plot((rec_data.time.rcc - 20e-3).*1e6,x_drive_current,'LineWidth',2);
        hold on
        plot((rec_data.time.rcc - 20e-3).*1e6,rec_drive_current,'LineWidth',2);

        xlim([0 300]);
        %ylim([0 1500]);
        set(gca,'fontsize',10);

        title('X-drive and Reconnection-drive Current');
        xlabel('Time [탎]');
        ylabel('Current [A]');

        hold on

        plot([ts_real ts_real], ylim, '-r','LineWidth',2);

        hold off
        drawnow;

       for frame_copy_index = 1:frame_rep
           frame_index = (ts_bdot-ts_bdot_begin)*frame_rep + frame_copy_index;
           frames{frame_index} = getframe(f);
       end

    end

    %Save video to file
    frames = cell2mat(frames);
    v = VideoWriter([langmuir_data_folder,'/profiles_video.mp4'],'MPEG-4');

    open(v);
    writeVideo(v,frames);
    close(v);
end