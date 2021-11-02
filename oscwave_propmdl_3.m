%% partial frequency (period) correction --> update freq each timepont & calculate predicted duration accordingly

clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;

% ======= make stim =======
stim = zeros(110,1)';
ioilist = [150, 120, 130, 140, 120];% [200, 200, 300, 600, 250];

for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 2;
    stim = [stim x];
end
fig = figure; 
set(fig,'Position', [ 1001         965        1085         374]);

subplot(3,1,1)
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);

% ======= init params =======
theta = 0;
freq_mdl = 5; % fs/samples << initial freq 
freq = freq_mdl;
mdl_ioi = fs/freq_mdl;

W = .6;
% WP = .5;

tbl = table();
a = 1;
k = 1;
t = Ts;

wave = [];
freqk = [];
thetak = [];
radphak = [];
radphak(1) = t*(2*pi*freq)+theta;
while 1     
    
    if stim(k) > 0  % if you see a stimulus (sound)

        if a == 1
            ioi = k; % ioi = n samples 
            prev_theta = theta;
        else
            ioi = k - tbl.k(a-1);
            prev_theta = thetak(k-1);
        end

        % -------------- update wave frequency  -------------- 
        freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
        
        % -------------- set required theta to match phase to prev --------------
        set_theta = prev_theta - t*2*pi*(freqnew - freq);
        
        % -------------- assign updated values  
        freq = freqnew;
        theta = set_theta; 
        
        contin_ok = input('continue? (1/0): ');
        if contin_ok ==0
            break;
        end
        
        % -------------- save interval info --------------
        tbl.ioi(a) = ioi;
        tbl.freq(a) = freq;
        tbl.k(a) = k;
        tbl.W(a) = W;
        a = a+1;

    end

    % -------------- save timepoint info --------------
    thetak(k) = theta;
    freqk(k) = freq;
    wave(k) = cos(t*(2*pi*freq)+theta); 
    
    % -------------- plot wave generated up to now --------------
    plot(wave, 'linewidth', 1.2, 'color', 'm');
    
    % increase time 
    t = t+Ts;
    if k == length(stim)
        break
    end    
    k = k+1;
end
EK_plotlabels('time', 'cos()', 'model wave',18);

subplot(3,1,2)
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(freqk, 'linewidth', 2, 'color', 'g'); 
EK_plotlabels('time', 'model frequency', '',18);

subplot(3,1,3)
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(thetak, 'linewidth', 2, 'color', 'g');
yline(0, '--k')
EK_plotlabels('time', 'model phase', '',18);



%% func

% function [phaneeded] = phase_finder (f2, pha2match, currpha, t)
% 
% 
% 
% radpha = t*2*pi*f2+currpha;
% 
% radphadif = pha2match-radpha;
% 
% 
% phaneeded = phaneeded/360;
% 
% end