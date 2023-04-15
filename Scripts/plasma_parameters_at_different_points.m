addpath('../Classes');
addpath('./Test');

%Plasma signals
plasma_source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_15_ECRH_Langmuir_Profile_10kHz_5V\oscilloscope\';
if ~exist('plasma_ST_profiles')
    load(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
%     plasma_ST_profiles.buildMDArrays;
%     plasma_ST_profiles.loadPpTimeserieCaFromFileList;
%     save(strcat(plasma_source_folder,'plasma_parameters_timeseries_list\','plasma_ST_profiles.mat'),'plasma_ST_profiles');
end

Nx = plasma_ST_profiles.Nx;
Ny = plasma_ST_profiles.Ny;
x_axis = plasma_ST_profiles.x_axis;
y_axis = plasma_ST_profiles.y_axis;

T_e_MDa = plasma_ST_profiles.T_e_MDa;
I_sat_i_MDa = plasma_ST_profiles.I_sat_i_MDa;
V_plasma_MDa = plasma_ST_profiles.V_plasma_MDa;

N = 3;
pp_timeseries_list = cell(N,1);
point_pos = nan(N,2);
%Center
point_pos(1,:) = [148,115];
%Ring
point_pos(2,:) = [70,70];
%Periphery
point_pos(3,:) = [20,60];

for n = 1:N
    if (~isnan(point_pos(n,1)) || ~isnan(point_pos(n,2)))
        pp_timeseries_list{n,1} = plasma_ST_profiles.getPpTimeseries(point_pos(n,1),point_pos(n,2),1);
    else
        fprintf('(x,y)-position is nan: no pp_timeseries selected\n');
    end
end

draw_T_e = true;
draw_n_e = true;
draw_V_plasma = true;

%Temperature profile
if draw_T_e
    %Average along time_axis (rise)
    T_e_rise = mode(T_e_MDa(:,:,1,:,1),4);
    %Average along time_axis (fall)
    T_e_fall = mode(T_e_MDa(:,:,1,:,2),4);
    %Average rise/fall
    T_e = (T_e_rise + T_e_fall)/2;
    %CSP smoothing
    csp_smoothing_param = 0.8;
    T_e = csp(T_e,csp_smoothing_param);
    %Isolines number
    n_isolines_T_e = 10;
    %Sanity check
    T_min = 0;
    T_max = 16;
    for nx = 1:Nx
        for ny = 1:Ny
            if (T_e(nx,ny) > T_max)
                T_e(nx,ny) = nan;
            end
        end
    end

    if(exist('f_T_e') && isvalid(f_T_e))
        if(exist('ax_T_e') && isvalid(ax_T_e))
            cla(ax_T_e);
        end
        clf(f_T_e);
    else
        f_T_e = figure;
    end
    figure(f_T_e);
    set(gcf,'Renderer','painters');
    ax_T_e = axes;
    [~,h_T_e_contourf] = contourf(ax_T_e,x_axis,y_axis,T_e',n_isolines_T_e);
    set(h_T_e_contourf,'Linecolor','k');
    caxis(ax_T_e, [T_min T_max]);
    set(ax_T_e,'fontsize',20);
    title(ax_T_e,'Electron temperature');
    xlabel(ax_T_e,'X position [mm]');
    ylabel(ax_T_e,'Y position [mm]');
    axis(ax_T_e,'equal');
    colormap(ax_T_e,'hot');
    h_T_e_cbar = colorbar(ax_T_e);
    ylabel(h_T_e_cbar, 'T_e [eV]');
    hold(ax_T_e,'on');
end

%Density profile
if draw_n_e
    %Average along time_axis (rise)
    I_sat_i_rise = mode(I_sat_i_MDa(:,:,1,:,1),4);
    %Average along time_axis (fall)
    I_sat_i_fall = mode(I_sat_i_MDa(:,:,1,:,2),4);
    %Average rise/fall
    I_sat_i = (I_sat_i_rise + I_sat_i_fall)/2;

    %USEFUL PHYSICAL CONSTANTS----------------------------------
    %Electron_charge [C]
    e = -1.6021766208 * 1e-19;
    %Boltzmann constant [J/K]
    k_B = 1.38064852 * 1e-23;
    %Atomic Mass Unit [kg]
    AMU = 1.660539040*1e-27;
    %-----------------------------------------------------------

    %Convert from I_sat_i to n_e
    %Probe tip area
    L_tip = 4.00*1e-3; %4.00 mm length
    r_tip = 0.25*1e-3; %0.25 mm radius
    probe_tip_area = L_tip*2*pi*r_tip;
    % Ion mass number
    ion_mass_number = 39.948;%Argon

    %Convert electron temperature(s) to SI units (Kelvin)
    T_e_SI = T_e.* ((-e)/k_B);
    %Convert ion mass from AMU to kg
    m_i = ion_mass_number*AMU;
    %Calculate electron charge density/ies [C*(m^-3)]
    rho_e = (I_sat_i./probe_tip_area) .* sqrt(2*pi*m_i./(k_B.*T_e_SI));
    %Calculate electron density [m^-3]
    n_e = rho_e./e;

    %CSP smoothing
    csp_smoothing_param = 0.8;
    n_e = csp(n_e,csp_smoothing_param);
    %Isolines number
    n_isolines_n_e = 10;
    %Sanity check
    n_min = 0;
    n_max = 5e18;
    for nx = 1:Nx
        for ny = 1:Ny
            if (I_sat_i(nx,ny) > n_max)
                I_sat_i(nx,ny) = nan;
            end
        end
    end

    if(exist('f_n_e') && isvalid(f_n_e))
        if(exist('ax_n_e') && isvalid(ax_n_e))
            cla(ax_n_e);
        end
        clf(f_n_e);
    else
        f_n_e = figure;
    end
    figure(f_n_e);
    set(gcf,'Renderer','painters');
    ax_n_e = axes;
    [~,h_n_e_contourf] = contourf(ax_n_e,x_axis,y_axis,n_e',n_isolines_n_e);
    set(h_n_e_contourf,'Linecolor','k');
    caxis(ax_n_e, [n_min n_max]);
    set(ax_n_e,'fontsize',20);
    title(ax_n_e,'Electron density');
    xlabel(ax_n_e,'X position [mm]');
    ylabel(ax_n_e,'Y position [mm]');
    axis(ax_n_e,'equal');
    colormap(ax_n_e,'parula');
    h_n_e_cbar = colorbar(ax_n_e);
    ylabel(h_n_e_cbar, 'n_e [m^{-3}]');
    hold(ax_n_e,'on');
end


%Draw points on profile
for n = 1:N
    hpT = plot(ax_T_e,point_pos(n,1),point_pos(n,2),'*b');
    hpT.LineWidth = 1;
    hpT.MarkerSize = 20;
    hpn = plot(ax_n_e,point_pos(n,1),point_pos(n,2),'*b');
    hpn.LineWidth = 1;
    hpn.MarkerSize = 20;
end

for n = 1:N
    time_axis_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.time_axis;
    V_signal_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.V_signal;
    V_ADC_step_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.V_ADC_step;
    I_signal_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.I_signal;
    I_smooth_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.I_smooth;
    I_smooth_exp_fit_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.I_smooth_exp_fit;
    I_CM_corrected_list{n} = pp_timeseries_list{n,1}.ivt_source.ivs_source_plasma.I_CM_corrected;
    time_axis_rise_list{n} = pp_timeseries_list{n,1}.ivt_source.total_time_axis_rise;
    time_axis_fall_list{n} = pp_timeseries_list{n,1}.ivt_source.total_time_axis_fall;
    NC_tot_list{n} = pp_timeseries_list{n,1}.ivt_source.NC_tot;
    nc = 10;
    V_axis_list{n} = pp_timeseries_list{n,1}.ivt_source.fall(nc).V_axis;
    I_of_V_list{n} = pp_timeseries_list{n,1}.ivt_source.fall(nc).I_of_V;
    I_of_V_fitted_list{n} = pp_timeseries_list{n,1}.I_of_V_fitted_ca{10,2};
    
    %V and I signals
    fig_ca{n,1} = figure;
    set(fig_ca{n,1}, 'NumberTitle', 'off','Name', ('V and I signals'));
    figure(fig_ca{n,1});
    
    set(fig_ca{n,1},'Renderer','painters');
    subplot(2,1,1)
    plot(time_axis_list{n},V_signal_list{n});
    subplot(2,1,2)
    plot(time_axis_list{n},I_CM_corrected_list{n});
    
    %V and I_CM_corrected signals
    fig_ca{n,2} = figure;
    set(fig_ca{n,2}, 'NumberTitle', 'off','Name', ('V and I_CM_corrected signals'));
    figure(fig_ca{n,2});
    set(fig_ca{n,2},'Renderer','painters');
    subplot(2,1,1)
    plot(time_axis_list{n},V_signal_list{n});
    subplot(2,1,2)
    plot(time_axis_list{n},I_CM_corrected_list{n});
    
    %IV-charateristic
    fig_ca{n,3} = figure;
    ta = time_axis_fall_list{n};
    time = ta(nc);
    set(fig_ca{n,3}, 'NumberTitle', 'off','Name', strcat('IV characteristic (with fit) (fall) at time = ',time*1000,' ms'));
    figure(fig_ca{n,3});
    set(fig_ca{n,3},'Renderer','painters');
    plot(V_axis_list{n},I_of_V_list{n},'*');
    hold
    plot(V_axis_list{n},I_of_V_fitted_list{n},'-');
    title(strcat('IV characteristic (with fit) (fall) at time = ',time*1000,' ms at point x = ',num2str(point_pos(n,1)),' mm y = ',num2str(point_pos(n,2)),' mm.'))
end





