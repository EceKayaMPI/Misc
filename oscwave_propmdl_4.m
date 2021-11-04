%% PROPMDL >> phase & freq corr OK. >> diff models with WP

clear all;clc;
close all;

fs = 1000;
Ts = 1/fs;

% ======= make stim =======
stim = zeros(100,1)';
% [200, 200, 200, 200, 200*6, 220]
% [200, 200, 300, 600, 250]
ioilist = [120, 110, 120, 200, 150];

for i = 1:length(ioilist)
    ons = zeros(ioilist(i),1)';
    ons(1) = 2;
    stim = [stim ons];
end
fig = figure;
set(fig,'Position', [1001 965 1085 374]);

% subplot(3,1,1)
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);



% W = .9;
% WP = 1;
weightlist = [0 .6 1];

for i = 1:3
% ======= init params =======
theta = 0;
mdl_freq = 6; % fs/samples << initial freq
freq = mdl_freq;
mdl_ioi = fs/mdl_freq; 

    W = 0;
    WP = weightlist(i);
    
    ons = 1;
    k = 1;
    t = Ts;
    
    wave = [];
    freqk = [];
    thetak = [];
    stimsamp = [];
    
    while 1
        
        if stim(k) > 0  % if you see a stimulus (sound)
            
            if ons == 1       % if it is the first one you saw
                % use your freq & phase values
                % because there isn't any interval yet. Just one sound.
                ioi = NaN;
                
            else % if it is the second or more , entrain...
                ioi = k - stimsamp(ons-1);
                prev_theta = thetak(k-1);
                
                
                % -------------- update wave frequency  (W) --------------
                stim_freq = fs/ioi;
                freqnew = (mdl_freq + (stim_freq - mdl_freq) * W); % avg freq
                %                 disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
                
                % -------------- set required theta to match phase to prev --------------
                set_theta = mod(prev_theta - t*2*pi*(freqnew - freq),2*pi);
                reset_theta = mod(2*pi*(1 - t*freqnew), 0);
                
                % -------------- update wave theta  (WP) --------------
                thetanew = (set_theta + (reset_theta - set_theta) * WP); % weighted theta correction
                
                % -------------- assign updated values
                freq = freqnew;
                theta = thetanew; % all theta is mod(theta)
                
                
            end % if first stim sound or not
            
            % -------------- save interval info --------------
            stimsamp(ons) = k;
            ons = ons+1;
        end
        
        % -------------- save timepoint info --------------
        thetak(k) = theta;
        freqk(k) = freq;
        wave(k) = cos(mod(t*(2*pi*freq),2*pi)+theta);%cos(t*(2*pi*freq)+theta);
        
        % increase time
        t = t+Ts;
        if k == length(stim)
            break
        end
        k = k+1;
    end
    
    subplot(3,1,i)
    plot(wave, 'linewidth', 1.2, 'color', 'm'); hold on;
    EK_xAxisMarker(find(stim>0), [0 0 1]);
    EK_plotlabels('time', 'wave', ['W-freq = ' num2str(W) ', W-theta = ' num2str(WP)], 18);
    
    
    
end

%% clean & change if stim or not

clear all;clc;
close all;

fs = 1000;
Ts = 1/fs;

% ======= make stim =======
stim = zeros(100,1)';
% [200, 200, 200, 200, 200*6, 220]
% [200, 200, 300, 600, 250]
ioilist = [120, 110, 120, 200, 150];

for i = 1:length(ioilist)
    ons = zeros(ioilist(i),1)';
    ons(1) = 2;
    stim = [stim ons];
end
fig = figure;
set(fig,'Position', [1001 965 1085 374]);



% W = .9;
% WP = 1;
weightlist = [0 .6 1];

for i = 1:3
    % ======= init params =======
    theta = 0;
    mdl_freq = 5; % fs/samples << initial freq
    freq = mdl_freq;
    mdl_ioi = fs/mdl_freq;
    
    W = 1;
    WP = weightlist(i);
    
    ons = 1;
    k = 1;
    t = Ts;
    
    wave = [];
    freqk = [];
    thetak = [];
    stimsamp = [];
    radvalk = [];
    
    while 1
        
        if stim(k) > 0  % if you see a stimulus (sound)
            
            if ons == 1       % if it is the first one you saw
                % use your freq & phase values
                % because there isn't any interval yet. Just one sound.
                ioi = NaN;
                
            else % if it is the second or more , entrain...
                ioi = k - stimsamp(ons-1);
                prev_theta = thetak(k-1);
                
                
                % -------------- update wave frequency  (W) --------------
                stim_freq = fs/ioi;
                freqnew = (mdl_freq + (stim_freq - mdl_freq) * W); % avg freq
                %                 disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
                
                % -------------- set required theta to match phase to prev --------------
                set_theta = prev_theta - t*2*pi*(freqnew - freq);
                %                 set_theta = radvalk(k-1) - t*2*pi*freq;
                reset_theta = 2*pi*(1 - t*freqnew);
                
                % -------------- update wave theta  (WP) --------------
                thetanew = (set_theta + (reset_theta - set_theta) * WP); % weighted theta correction
                
                % -------------- assign updated values
                freq = freqnew;
                theta = thetanew;
                
                
            end % if first stim sound or not
            
            % -------------- save interval info --------------
            stimsamp(ons) = k;
            ons = ons+1;
        end
        
        % -------------- save timepoint info --------------
        thetak(k) = theta;
        freqk(k) = freq;
        wave(k) = cos(t*(2*pi*freq)+theta);
        radvalk(k) = t*(2*pi*freq)+theta;
        
        % increase time
        t = t+Ts;
        if k == length(stim)
            break
        end
        k = k+1;
    end
    
    subplot(3,1,i)
    plot(wave, 'linewidth', 1.2, 'color', 'm'); hold on;
    EK_xAxisMarker(find(stim>0), [0 0 1]);
    EK_plotlabels('time', 'wave', ['W-freq = ' num2str(W) ', W-theta = ' num2str(WP)], 18);
    
    
    
end
%
% subplot(3,1,1)
% plot(wave, 'linewidth', 1.2, 'color', 'm'); hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% EK_plotlabels('time', 'wave', ['period correction weight = ' num2str(W) ', phase correction weight = ' num2str(WP)], 18);
%
%
% subplot(3,1,2)
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% plot(freqk, 'linewidth', 2, 'color', 'g');
% EK_plotlabels('time', 'model frequency', '',18);
%
% subplot(3,1,3)
% hold on;
% EK_xAxisMarker(find(stim>0), [0 0 1]);
% plot(radvalk, 'linewidth', 2, 'color', 'g');
% yline(0, '--k')
% EK_plotlabels('time', 'model radval', '',18);


%% NO phase reset, but entr decay (W decays)

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
mdl_ioi = fs/mdl_freq;

W = 0;
W_init = 1;
% WP = .5;

tbl = table();
ons = 1;
k = 1;
t = Ts;

wave = [];
freqk = [];
thetak = []; 
W_time = [];
W_decaycons = .01;
onslistk = [];

while 1
    
 
    if stim(k) > 0  % if you see a stimulus (sound)
       
        % check onsets & find stim frequency
        if ons > 1 % if it is first onset, do nothing
            
            ioi = k - onslistk(ons-1);
            
            % maximize entrainment weight
            W = W_init;
        
        end
  
        % save interval & onset info
        onslistk(ons) = k;
        ons = ons+1;   
    end    

    if k == 1
        prev_theta = theta;
    else
        prev_theta = thetak(k-1);
    end
    
    if W > 0
        % -------------- update wave frequency  --------------
        freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        

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
    wave(k) = cos(mod(t*(2*pi*freq),2*pi)+mod(theta,2*pi));% cos(t*(2*pi*freq)+theta);
    
    % -------------- plot wave generated up to now --------------
%     plot(wave, 'linewidth', 1.2, 'color', 'm');
    
    % increase time
    t = t+Ts;
    if k == length(stim)
        break
    end
    k = k+1;
end

fig = figure;
set(fig,'Position', [ 1001         965        1085         374]);
sgtitle({['Model freq = ' num2str(mdl_freq)],...
          ['W exponential decay (W = W - W * ' num2str(W_decaycons) ')']  }, 'fontsize', 18);


subplot(3,1,1) % stim 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(wave, 'linewidth', 2, 'color', 'm');
EK_plotlabels('time', 'cos()', 'wave',18);

subplot(3,1,2)  % freqs 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(freqk, 'linewidth', 2, 'color', 'g');
EK_plotlabels('time', 'frequency', 'model frequency',18);

subplot(3,1,3) % W 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(W_time, 'linewidth', 2, 'color', 'g');
yline(0, '--k')
EK_plotlabels('time', 'W', 'freq correction weight',18);

%% entr decay ++

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
freqk = [];
thetak = []; 
W_time = [];
W_decaycons = .01;
onslistk = [];

while 1
    
 
    if stim(k) > 0  % if you see a stimulus (sound)
       
        % check onsets & find stim frequency
        if ons > 1 % if it is first onset, do nothing        
            ioi = k - onslistk(ons-1);
            % maximize entrainment weight
            W = W_init;
        
        end
  
        % save interval & onset info
        onslistk(ons) = k;
        ons = ons+1;   
    end    

    if k == 1 % first sample 
        prev_theta = theta;
    else
        prev_theta = thetak(k-1);
    end
    
    if W > 0
        % -------------- update wave frequency  --------------
%         freqnew = fs/(mdl_ioi + (ioi - mdl_ioi) * W); % avg ioi
        stim_freq = fs/ioi;
        
   
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
    wave(k) = cos(mod(t*(2*pi*freq),2*pi)+mod(theta,2*pi));% cos(t*(2*pi*freq)+theta);
    
    % -------------- plot wave generated up to now --------------
%     plot(wave, 'linewidth', 1.2, 'color', 'm');
    
    % increase time
    t = t+Ts;
    if k == length(stim)
        break
    end
    k = k+1;
end

fig = figure;
set(fig,'Position', [ 1001         965        1085         374]);
sgtitle({['Model freq = ' num2str(mdl_freq)],...
          ['W exponential decay (W = W - W * ' num2str(W_decaycons) ')']  }, 'fontsize', 18);


subplot(3,1,1) % stim 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(wave, 'linewidth', 2, 'color', 'm');
EK_plotlabels('time', 'cos()', 'wave',18);

subplot(3,1,2)  % freqs 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(freqk, 'linewidth', 2, 'color', 'g');
EK_plotlabels('time', 'frequency', 'model frequency',18);

subplot(3,1,3) % W 
hold on;
EK_xAxisMarker(find(stim>0), [0 0 1]);
plot(W_time, 'linewidth', 2, 'color', 'g');
yline(0, '--k')
EK_plotlabels('time', 'W', 'freq correction weight',18);
