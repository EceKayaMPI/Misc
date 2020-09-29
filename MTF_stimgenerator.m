
clear,clc;

N = 34;
from = 1500;
to = 100;
stim_IOIs = [50 linspace(from,to,N)];

Fs = 44100;
stim = audioread('Woodblock3.wav'); % generate sound stimuli
stim_dur_ms = 100; % in milisecs
stim = EK_cut_and_fade (stim,Fs,stim_dur_ms);
stim = stim.* .2;
stim_dur_sec = length(stim)/Fs;

stim_seq = [];

for n = 1:N

IOI_ms = stim_IOIs(n);
    
ISI = zeros(round(Fs*(IOI_ms/1000 - stim_dur_sec)),2).';

stim_seq = [stim_seq stim ISI];

end

filename = sprintf('MTF_%d_to_%dms_%dtaps.wav', from, to, N);
audiowrite (filename, stim_seq(1,:), Fs);