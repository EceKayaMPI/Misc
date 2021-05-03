

stim_ms_list = 400:100:800;
stimnames = {};
HWE_INFO = struct(); 
HWE_INFO.stimlist = [];

[stimsound, Fs] = audioread('Woodblock3.wav');
[stim] = cut_and_fade (stimsound, Fs, 60, []);
% figure;plot(stim')

for i = 1:length(stim_ms_list)
IOI_ms = stim_ms_list(i);
[stim_seq, Fs] = isochron_seq (10, stim, Fs, IOI_ms, 7);

wavname = sprintf('%d_ms.wav',IOI_ms);
% audiowrite(wavname, stim_seq', Fs);
stimnames{i} = wavname;
end

k = 1;
for ii = 1:length(stim_ms_list)
    
    for jj = 1:length(stim_ms_list)
        
        if ii~=jj
            HWE_INFO.stimlist(k).stim_1 = stimnames{ii};
            HWE_INFO.stimlist(k).stim_2 = stimnames{jj};
            k = k+1;
        end
    end
end

filename = '/Users/ece.kaya/JS/JATOS/study_assets_root/HWE/HWE_PID.json';
savejson('',HWE_INFO,filename);
%% fun

function [stim_seq, Fs] = isochron_seq (sildur_ms, stim, Fs, IOI_ms, nStims)

    if IOI_ms < 300
        nEmp = nEmp+1;
    end
    
    silenceXms = zeros(round(Fs*(sildur_ms/1000)),2).';

    stim_dur_sec = length(stim)/Fs;
    ISI = zeros(round(Fs*(IOI_ms/1000 - stim_dur_sec)),2).'; 
	empIOI = zeros(round(Fs*(IOI_ms/1000)),2).';
    
    stim_seq = [stim ISI];
    stim_seq = repmat(stim_seq,1,nStims);
    
    stim_seq = [silenceXms stim_seq];
    
end

function [newstim] = cut_and_fade (sound, Fs, stim_dur_ms, fadein_ms)

stim_dur_sec = stim_dur_ms/1000; % in seconds >> because turning it into a sound wave will require so.
stim = sound(1:round(stim_dur_sec*Fs),1:2).'; % don't forget to put .' (or use transpose) for sound stuff.
fade_samp = round(stim_dur_sec*Fs);
fade_frame = linspace(1,0,fade_samp);
newstim = stim.*fade_frame;

% fadein_samp = round(fadein_ms/1000*Fs);
% fadein_frame = linspace(0,1,fadein_samp);
% newstim(:,1:fadein_samp) = newstim(:,1:fadein_samp).*fadein_frame;
end