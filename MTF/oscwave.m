
figure;

dur = 3;
f = 6;
fs = 44100;  
Ts = 1/44100;              % sampling period
% ncycs = f/dur;          % number of cycles

ncycs = 10;             %have to give int

% dur = f / ncycs;
% Tvec = 1/fs:Ts:dur;     % time vector. <<< changed to ones(..)

%--------single freq
test = [];
fi  = f;
for ii = 1:ncycs
    cycdur = 1/fi;
    cycT   = Ts:Ts:cycdur;
    test = [test cos(2*pi*fi*cycT+pi)];
end

% subplot(2,1,1)
plot(test); hold on;

testm = [];
flist = randi(2*f,ncycs,1);
for ii = 1:ncycs
    cycdur = 1/flist(ii);
    cycT   = Ts:Ts:cycdur;
    testm = [testm cos(2*pi*flist(ii)*cycT+pi)];
end

% subplot(2,1,2)
plot(testm)

%% 
clear; clc; close all;

fs = 44100;  
Ts = 1/fs;  
sinwa = [];

stim_ioi = 500;
gap_int = 3;

stim_seq = [repmat(stim_ioi,5,1)' stim_ioi*gap_int stim_ioi];
% stim_seq = [repmat(500,5,1)' 1500 500];

Wpha = .5;
Wper =  .5;

% pha = 0;            % phase corr term
% per = 500;          % period corr term
% C = 0;              % temporal contrast
pha = NaN(length(stim_seq),1)';
pha(1) = 0;

per = NaN(length(stim_seq),1)';
per(1) = 500;



chist = [];

% model stim seq loop -------------------
for i = 1:5%length(stim_seq)-1



    pha(i+1) = (1 - Wpha)*C;
    per(i+1) = (1 + Wper*C)*per(i);

    chist(i) = C;
    
end





stim_seq_times = [0 cumsum(stim_seq)];
per_seq_times = [0 cumsum(per)];

figure;
plot(stim_seq_times, ones(length(stim_seq_times),1) , '.k', 'MarkerSize', 30)
hold on;
% plot(sinwa, 'linewidth', 1.5)

plot(per_seq_times, ones(length(per_seq_times),1)+.5 , '.r', 'MarkerSize', 30)

xticks(stim_seq_times); grid on; 
yticks([1 1.5]); yticklabels({'stim', 'model'});
set(gca, 'YGrid', 'off', 'xlim', [-100 max(stim_seq_times)+100]);
EK_plotlabels('timepoints', '', ['C_{final} = ' num2str(C)], 15);
set(gcf,'Position',[1500 500 500 100])

