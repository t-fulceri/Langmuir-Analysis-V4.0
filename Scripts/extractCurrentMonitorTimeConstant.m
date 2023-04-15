%Thi script extracts the time constant tau which models the current monitor
%"capacitor charging"-like effect

%Tiziano Fulceri 2019-02-21

%Optimized to work with langmuir measurements with the following parameters:
% f = 10 kHz
% V_max = 200 V
% time interval = 10 ms
% Probe position = (148, 115) (Langmuir probe at the center)
% ECRH plasma ON

clear all
close all

addpath('../');

%Read the plasma measurement
source_folder = 'C:\Users\tzf\Desktop\VINETA II Lab\tzf_local\2019_01_16_ECRH_Langmuir_CentralPosition_10kHz_10V_DeamplifiedBy200\'
settings = struct;
settings.voltage_sign = 1;
settings.current_sign = -1;
settings.voltage_deamplification = 200;

%Number of repetition to be averaged
N = 10;
%Array of measurements
filename_arr = cell(N,1);
ivs_arr = cell(N,1);
for n = 1:N
    filename_arr{n,1} = strcat('meas_',num2str(n-1,'%04.f'),'.h5');
    ivs_arr{n,1} = IvSignals;
    ivs_arr{n,1} = ivs_arr{n,1}.assignDAQ_settings(settings);
    ivs_arr{n,1}.assignPropertiesFromHD5File(source_folder,filename_arr{n,1});
end

smooth_param = 1000;

time_axis = ivs_arr{1,1}.time_axis;
length = ivs_arr{1,1}.signal_length;
I_arr = zeros(length,N);
I_smooth_arr = zeros(length,N);
for n = 1:N
    I_arr(:,n) = ivs_arr{n,1}.I_signal;
    I_smooth_arr(:,n) = smooth(I_arr(:,n),smooth_param,'sgolay');
end

%Take the maxima of the smoothed values: they will serve as the amplitude
%input parameter for the fitting procedure
I_smooth_max_arr = max(I_smooth_arr);

%Time constant input parameter (based on observations of current signal)
tau_input = 1e-3;

fitresult_arr = cell(N,1);
gof_arr = cell(N,1);

for n = 1:N
    [xData, yData] = prepareCurveData( time_axis, I_smooth_arr(:,n) );

    % Set up fittype and options.
    ft = fittype( 'A*exp(-t/tau)', 'independent', 't', 'dependent', 'I' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [I_smooth_max_arr(n) tau_input];

    % Fit model to data.
    [fitresult_arr{n,1}, gof_arr{n,1}] = fit( xData, yData, ft, opts );
end

tau_arr = zeros(N,1);

for n = 1:N
    tau_arr(n) = fitresult_arr{n,1}.tau;
end

tau = mean(tau_arr)
tau_std = std(tau_arr)





