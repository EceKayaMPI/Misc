
clear all;clc;
close all;

fs = 1000;
Ts = 1/fs;

% ======= make stim =======
stim = zeros(200,1)';
ioilist = [100, 100, 100, 100, 100*6,100];% [200, 200, 300, 600, 250];

for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 2;
    stim = [stim x];
end


% ======= init params =======
theta = 0;
mdl_freq = 4; % fs/samples << initial freq
freq = mdl_freq;
% mdl_ioi = fs/mdl_freq;

W = 0;
W_init = 1;
% WP = .5;

tbl = table();
ons = 1;
k = 1;
t = Ts;

wave = [];
phasek = [];
freqk = [];
thetak = []; 
W_time = [];
W_decaycons = .01;


while 1
    
 
    if stim(k) > 0  % if you see a stimulus (sound)
       disp(phasek(k-1));
       
       corr_freq_percent = ( 360 - rad2deg(phasek(k-1)) )./360;
       % .5-.5 >> make faster
       break;

    end    
 
    if W > 0
        % -------------- update wave frequency  --------------
%         freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        corr_freq = mdl_freq + mdl_freq*corr_freq_percent;
        
   
        freqnew = (mdl_freq + (stim_freq - mdl_freq) * W);

        % -------------- set required theta to match phase to prev --------------
        set_theta = prev_theta - t*2*pi*(freqnew - freq);

        % -------------- assign updated values
        freq = freqnew;
        theta = set_theta;

        % -------------- entr decay --------------
        W = W - W * W_decaycons;
    end
    
    % -------------- save timepoint info --------------
    W_time(k) = W;
    thetak(k) = theta;
    freqk(k) = freq;
    wave(k) = cos(t*(2*pi*freq)+theta);
    phasek(k) = t*(2*pi*freq)+theta;
    
    % -------------- plot wave generated up to now --------------
%     plot(wave, 'linewidth', 1.2, 'color', 'm');
    
    % increase time
    t = t+Ts;
    if k == length(stim)
        break
    end
    k = k+1;
end

% fig = figure;
% set(fig,'Position', [ 1001         965        1085         374]);
% sgtitle({['Model freq = ' num2str(mdl_freq)],...
%           ['W exponential decay (W = W - W * ' num2str(W_decaycons) ')']  }, 'fontsize', 18);
% 
% 
% subplot(3,1,1) % stim 
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% plot(wave, 'linewidth', 2, 'color', 'm');
% EK_plotlabels('time', 'cos()', 'wave',18);
% 
% subplot(3,1,2)  % freqs 
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% plot(freqk, 'linewidth', 2, 'color', 'g');
% EK_plotlabels('time', 'frequency', 'model frequency',18);
% 
% subplot(3,1,3) % W 
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% plot(W_time, 'linewidth', 2, 'color', 'g');
% yline(0, '--k')
% EK_plotlabels('time', 'W', 'freq correction weight',18);
