% figure
% set(gcf,'Renderer','painters');
% movegui(gcf,'northwest');
% subplot(2,1,1);
% plot(time_axis,V_signal_plasma);
% subplot(2,1,2);
% plot(time_axis,I_signal_plasma);
% 
% figure
% set(gcf,'Renderer','painters');
% movegui(gcf,'north');
% subplot(2,1,1);
% plot(time_axis,V_signal_vacuum);
% subplot(2,1,2);
% plot(time_axis,I_signal_vacuum);
% 
% figure
% set(gcf,'Renderer','painters');
% movegui(gcf,'west');
% subplot(2,1,1);
% plot(time_axis,V_signal_plasma);
% subplot(2,1,2);
% plot(time_axis,I_signal_plasma_CM_corrected);

% figure
% set(gcf,'Renderer','painters');
% movegui(gcf,'center');
% subplot(2,1,1);
% plot(time_axis,V_signal_plasma);
% subplot(2,1,2);
% plot(time_axis,I_signal_plasma_CM_corrected - I_signal_vacuum);

% figure
% set(gcf,'Renderer','painters');
% movegui(gcf,'east');
% plot(time_axis,I_signal_plasma - I_signal_vacuum);
% hold on
% plot(time_axis,smooth(I_signal_plasma - I_signal_vacuum,1000,'lowess'));
% plot(time_axis,I_smooth_exp_fit);

