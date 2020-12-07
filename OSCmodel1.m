clear,clc;
% parameters
Ts = 1/1000;
sig = [];


Wphi = .5;              % 0 = no phase correction, 1 = full phase reset 
Wp = .4;                % 0 = no period corr. 1 = complete period corr where P matches curr IOI 

pha = 0;            % phase corr term
per = 600;            % period corr term
C = 0;                    % temporal contrast

IOI = 600;                  % IOI of the stimulus
dev_perc = .1;

stim_seq = [repmat(IOI,1,7) 540 1200 486]; 
on_time = repmat(IOI,1,11);

C_all = [];

for i = 1:length(stim_seq)

if mod((pha(i) + stim_seq(i)/per(i)),1) > .5
    C = mod((pha(i) + stim_seq(i)/per(i)),1) - 1;
else
    C = mod((pha(i) + stim_seq(i)/per(i)),1);
end

pha(i+1) = (1 - Wphi)*C;
per(i+1) = (1 + Wp*C)*per(i);

C_all(i) = C;


    cycdur = per(i);
    fi = 1/cycdur;
    cycT   = Ts:Ts:cycdur;
    wavez{i} = cos(2*pi*fi*cycT+pha(i));
    sig = [sig wavez{i}];  

end


disp(C);





figure;
scatter(cumsum(stim_seq), ones(1,length(cumsum(stim_seq))), 'filled', 'k'); hold on;
scatter(cumsum(on_time), ones(1,length(cumsum(on_time))), 'k');


figure;
scatter(cumsum(stim_seq).*1000, ones(1,length(cumsum(stim_seq))), 'filled', 'k'); hold on;
plot(sig, 'LineWidth', 2)
ylim ([-1.1 1.1])
