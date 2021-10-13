% 
close all; clear; clc;

kappa = .5;
t_info = table();

figure;
subplot(2,2,1)
k = 1;
for t = 1:.1:2*pi % sampling rate  
    
    theta = t;
    theta_fast = 1.3*t;
    
    t_info.t(k) = t;
    t_info.theta(k) = theta;
    t_info.theta_fast(k) = theta_fast;
    k = k+1;

    make_dots(theta); hold on;
    make_dots(theta_fast); hold off;
    drawnow
end

subplot(2,2,2)
plot(t_info.t, sin(t_info.theta)); hold on;
plot(t_info.t, sin(t_info.theta_fast)); 



%%
% close all; clear; clc;
% 
% timestep = .5;  % step that timer func evaluates (sampling rate)
% timeInc = .1;   % increments in samples at every function call 
% 
% figure;
% 
% for t = 1:.1:10
%     
% %     theta = mod(t, 2*pi)+pi;
% %     theta_fastosc = mod(1.3*t, 2*pi)+pi;
%     
%     theta = t;
%     theta_fastosc = 1.3*t;
% 
%     make_dots(theta); hold on;
%     make_dots(theta_fastosc); hold off;
%     drawnow
% end
            
%% functions

% function kuramoto(thetaVec, k, 

function make_dots(theta)
    plot(cos(theta), sin(theta), "Marker",".","MarkerSize",30); 
    set(gca,'xlim', [-1 1], 'ylim', [-1 1]);
end


