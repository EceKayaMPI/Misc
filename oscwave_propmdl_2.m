
%% full frequency (period) correction 
clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;

wave = [];

ioilist = [100, 100, 100, 100, 100];
stim = zeros(200,1)';
for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 1;
    stim = [stim x];
end

figure; stem(stim); hold on;
% make stim VECTOR 


k = 1;
t = Ts;

freq = 5; % fs/samples << initial freq 

tbl = table();
a = 1;

while 1     
    
    if stim(k) == 1 % if you see a stimulus (sound)

        if a == 1
            ioi = k;
        else
            ioi = k - tbl.k(a-1);
        end

        
        freqnew = fs/ioi;
        disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
        freq = freqnew;        
        
        contin_ok = input('continue? (1/0): ');
        if contin_ok ==0
            break;
        end
        
        tbl.ioi(a) = ioi;
        tbl.freq(a) = freq;
        tbl.k(a) = k;
        a = a+1;
    end
    
    
    
    % calculates by adding wave points, calc by multiplying time value
    % (t)
    wave(k) = cos(t*(2*pi*freq)+0);
    plot(wave);
    
    
    t = t+Ts;
    if k == length(stim)
        break
    end
    
    k = k+1;
end

%% partial frequency (period) correction 

clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;

wave = [];

ioilist = [100, 120, 130, 140, 120];
stim = zeros(200,1)';
for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 1;
    stim = [stim x];
end

figure; bar(stim); hold on;
% make stim VECTOR 


k = 1;
t = Ts;

% freq = 5; % fs/samples << initial freq 
mdl_ioi = 220;
freq = fs/mdl_ioi;

tbl = table();
a = 1;

W = .999;
WP = .5;

pha = pi/2;

% %         % McAuley model equations
% %         pha(i+1) = (1 - Wpha)*C;
% %         per(i+1) = (1 + Wper*C)*per(i);

while 1     
    
    if stim(k) == 1 % if you see a stimulus (sound)

        if a == 1
            ioi = k;
%             pha = 0;
        else
            ioi = k - tbl.k(a-1);

        end

        % -------------- update wave frequency  -------------- 
        freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
        
        % [pha3] = EK_phase_matcher (prev_freq, curr_freq, prev_phase, curr_phase, t)
        pha = EK_phase_matcher (prev_freq, curr_freq, prev_phase, curr_phase, t)
        
        
        freq = freqnew;        
        
        contin_ok = input('continue? (1/0): ');
        if contin_ok ==0
            break;
        end
        

        
        tbl.pha(a) = pha;
        tbl.ioi(a) = ioi;
        tbl.freq(a) = freq;
        tbl.k(a) = k;
        tbl.W(a) = W;
        a = a+1;
    end
    
    
    
    % calculates by adding wave points, calc by multiplying time value
    % (t)
    wave(k) = cos(t*(2*pi*freq)+pha);
    plot(wave);
    
    
    t = t+Ts;
    if k == length(stim)
        break
    end
    
    k = k+1;
end

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

subplot(2,1,1)
plot(stim -1, 'linewidth', 2); hold on;

% ======= init params =======
pha = 0;
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
while 1     
    
    if stim(k) > 0  % if you see a stimulus (sound)

        if a == 1
            ioi = k; % ioi = n samples 
            prev_pha = pha;
        else
            ioi = k - tbl.k(a-1);
            prev_pha = tbl.pha(a-1);
        end

        % -------------- update wave frequency  -------------- 
        freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
        
        
        prev_pha = phasek(k-1);
        set_pha = 
%         wave(k) = cos(t*(2*pi*freq)+pha); 
        
        freq = freqnew;
        pha = set_pha; 
        
        contin_ok = input('continue? (1/0): ');
        if contin_ok ==0
            break;
        end
        
        
        
        
        % -------------- save interval info --------------
        tbl.pha(a) = pha;
        tbl.ioi(a) = ioi;
        tbl.freq(a) = freq;
        tbl.k(a) = k;
        tbl.W(a) = W;
        a = a+1;

                
        
    end

    % -------------- save timepoint info --------------
    % calculates by adding wave points, calc by multiplying time value
    % (t)
    freqk(k) = freq;
    % this is radian: t*(2*pi*freq)
    % we add +pha which shifts it
    % take inside cos
    phasek(k) = t*(2*pi*freq)+pha;
    
    % here > replace + pha with 'phase correction'
    % phase corr: calculate diff in phase (what should have been 

    wave(k) = cos(t*(2*pi*freq)+pha); 
    
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

subplot(2,1,2)
plot(stim/2*(max(freqk)+1),'linewidth', 2); hold on;
plot(freqk, 'linewidth', 2);
EK_plotlabels('time', 'model frequency', '',18);

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