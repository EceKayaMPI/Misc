
clear all;clc;
close all;

fs = 1000;
Ts = 1/fs;

% ======= make stim =======
stim = zeros(200,1)';
ioilist = [100, 100, 100, 100, 100*6];% [200, 200, 300, 600, 250];

for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 2;
    stim = [stim x];
end

% ======= init params =======
theta = 0;
mdl_freq = 4; % fs/samples << initial freq
freq = mdl_freq;

k = 1;
t = Ts;

wave = [];
phasek = [];

while 1
    
    
    if stim(k) > 0  % if you see a stimulus (sound)  
        
        % check your phase
        currpha = (t*(2*pi*freq)+theta);
        
        
        
        
       
        
    end
    
    
    % -------------- save timepoint info --------------
    wave(k) = cos(t*(2*pi*freq)+theta);
    phasek(k) = t*(2*pi*freq)+theta;

    % increase time
    t = t+Ts;
    if k == length(stim)
        break
    end
    k = k+1;
    
end

figure;
subplot(2,1,1)
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(wave, 'linewidth', 2, 'color', 'm');
EK_plotlabels('time', 'cos()', 'wave',18);

subplot(2,1,2)
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(phasek, 'linewidth', 2, 'color', 'm');
EK_plotlabels('time', 'phase', 'phase',18);

