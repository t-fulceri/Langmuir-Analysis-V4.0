%Plot the v_ALFVEN as a function of background magnetic field and ion density
%Plot the v_e,thermal as a function of electron temperature
%Tiziano Fulceri 2019-05-29


%ALFVEN VELOCITY
%Gas mass number
%Argon
A_ion = 39.948;

B_min = 0e-3;
B_max = 50e-3;
B_arr_length = 200;
B_arr = linspace(B_min,B_max,B_arr_length);

n_min = 1e17;
n_max = 1e18;
n_arr_length = 100;
n_arr = linspace(n_min,n_max,n_arr_length);

B_mat = repmat(B_arr',[1,n_arr_length]);
n_mat = repmat(n_arr,[B_arr_length,1]);

v_A = zeros(B_arr_length,n_arr_length);
v_A(:,:) = PlasmaPhysics.alfven_velocity(B_mat(:,1),A_ion,n_arr(1,:),false);


%ELECTRON THERMAL SPEED
T_e_min = 0;
T_e_max = 20;
T_e_arr_length = 100;
T_e_arr = linspace(T_e_min,T_e_max,T_e_arr_length);

close all

set(0, 'DefaultFigureRenderer', 'opengl');

figure
[c,h] = contourf(1e3.*B_arr,1e-16.*n_arr,1e-5.*v_A',30);
clabel(c,h);
movegui(gcf,'northwest');
title('Alfven velocity as a function of background magnetic field and ion density')
xlabel('B_0 [mT]');
ylabel('n_i[10^{16} m^{-3}]');
h_cbar = colorbar;
ylabel(h_cbar, 'v_A(B_0,n_i) [10^{5} m/s]');

figure
plot(T_e_arr,1e-5.*PlasmaPhysics.electron_thermal_speed(T_e_arr,false));
movegui(gcf,'north');
title('Electron thermal speed as a function of electron temperature')
xlabel('T_e [eV]');
ylabel('v_{e_th} [10^{5} m/s]');
