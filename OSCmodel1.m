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

end_early = [repmat(IOI,1,7) (IOI-IOI*dev_perc) IOI*2 (IOI-IOI*dev_perc)]; 
end_late = [repmat(IOI,1,7) (IOI+IOI*dev_perc) IOI*2 (IOI+IOI*dev_perc)];
beg_early = [repmat(IOI,1,7) IOI*2-IOI*dev_perc (IOI+IOI*dev_perc)];
beg_late = [repmat(IOI,1,7) IOI*2+IOI*dev_perc (IOI-IOI*dev_perc)];
on_time = [repmat(IOI,1,7) IOI*2 IOI];

stim_seq = beg_late;

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

%     fi  = ms2Hz(per(i));
%     cycdur = 1/fi;
%     cycT   = Ts:Ts:cycdur;
%     wavez{i} = cos(2*pi*fi*cycT+pha(i));
%     sig = [sig cos(2*pi*fi*cycT+pha(i))];  
%     plot(sig)


    cycdur = per(i);
    fi = 1/cycdur;
    cycT   = Ts:Ts:cycdur;
    wavez{i} = cos(2*pi*fi*cycT+pha(i));
    sig = [sig wavez{i}];  
    plot(sig)


end


disp(C);
% Luce_condition = exp(gamma*C_condition)/sum(exp(gamma*C_all_conditions))







% 
% figure;
% subplot(2,1,1)
% scatter(cumsum(stim_seq), ones(1,length(cumsum(stim_seq))), 'filled', 'k'); hold on;
% scatter(cumsum(on_time), ones(1,length(cumsum(on_time))), 'k');
% 
% subplot(2,1,2)
% plot(test)


%%
    fi  = ms2Hz(per(i));
    cycdur = 1/fi;
    cycT   = Ts:Ts:cycdur;

% x = cos(2*pi*fi*cycT+pha(i));
% plot(x)

%% func

function EK_dotplot(x)
scatter(x, ones(1,length(x)),'LineWidth',1);
set(gca,'LineWidth',0.75);
end