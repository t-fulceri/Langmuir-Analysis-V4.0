V_tf_final = V_axis(1:knee_ind);
I_tf_final = I_of_V(1:knee_ind);

[xData, yData] = prepareCurveData( V_tf_final, I_tf_final );

% Set up fittype and options.
ft = fittype( 'I_sat_i + k*(V - V_plasma)+ I_sat_e*exp(b*(V-V_plasma))', 'independent', 'V', 'dependent', 'I' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -Inf -Inf 0 0];
opts.Robust = 'Bisquare';
opts.StartPoint = [I_sat_e_estimate I_sat_i_estimate V_axis(knee_ind) inv_T_e_estimate sheath_coeff];
opts.Upper = [max(I_sat_e_estimate_bounds) 0  max(V_plasma_estimate_bounds) max(inv_T_e_estimate_bounds) max(sheath_coeff_bounds)];
opts.Lower = [min(I_sat_e_estimate_bounds) min(I_sat_i_estimate_bounds) min(V_plasma_estimate_bounds) min(inv_T_e_estimate_bounds) min(sheath_coeff_bounds)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts )

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'I_tf_final vs. V_tf_final', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel V_tf_final
ylabel I_tf_final
grid on