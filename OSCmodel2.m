
clear,clc;
% parameters

Ts = 1/1000;
sig = [];

Wphi = .5;              % 0 = no phase correction, 1 = full phase reset
Wp = .4;                % 0 = no period corr. 1 = complete period corr where P matches curr IOI

pha = 0;            % phase corr term
per = 600;            % period corr term
C = 0;                    % temporal contrast

IOIlist = [540, 600, 660];

appdx = struct();

for ioi = 1:1
%     IOI = IOIlist(ioi);
    IOI = 600;
    dev_perc = .1;
    
    end_early = [repmat(IOI,1,7) (IOI-IOI*dev_perc) IOI*2 (IOI-IOI*dev_perc)];
    end_late = [repmat(IOI,1,7) (IOI+IOI*dev_perc) IOI*2 (IOI+IOI*dev_perc)];
    beg_early = [repmat(IOI,1,8) IOI*2-IOI*dev_perc (IOI+IOI*dev_perc)];
    beg_late = [repmat(IOI,1,8) IOI*2+IOI*dev_perc (IOI-IOI*dev_perc)];
    on_time = [repmat(IOI,1,8) IOI*2 IOI];

    all_condit_seqs = vertcat(end_early, on_time, end_late, beg_early, on_time, beg_late);

    appdx(ioi).stimz = all_condit_seqs;

    for condit = 1:6
    
    stim_seq = all_condit_seqs(condit,:);
    
    % model stim seq loop -------------------
    for i = 1:length(stim_seq)
        
        if mod((pha(i) + stim_seq(i)/per(i)),1) > .5
            C = mod((pha(i) + stim_seq(i)/per(i)),1) - 1;
        else
            C = mod((pha(i) + stim_seq(i)/per(i)),1);
        end
        
        pha(i+1) = (1 - Wphi)*C;
        per(i+1) = (1 + Wp*C)*per(i);
        
    end
    
    appdx(ioi).C_fin(condit) = C;
    end
    
end % IOIlist(k) loop end



% figure;
% subplot(2,1,1)
% scatter(cumsum(stim_seq), ones(1,length(cumsum(stim_seq))), 'filled', 'k'); hold on;
% scatter(cumsum(on_time), ones(1,length(cumsum(on_time))), 'k');
% 
% subplot(2,1,2)
% plot(sig)
% 
% 
% figure
% scatter(cumsum(stim_seq).*1000, ones(1,length(cumsum(stim_seq))), 'filled', 'k'); hold on;
% plot(sig, 'LineWidth', 2)
% ylim ([-1.1 1.1])
%% func

function EK_dotplot(x)
scatter(x, ones(1,length(x)),'LineWidth',1);
set(gca,'LineWidth',0.75);
end