
%% oscillator freezes during intervals that double its period

clear; clc; close all;

% fs = 44100;
% Ts = 1/fs;
% sinwa = [];

stim_ioi = 500;
stim_comparison = stim_ioi + stim_ioi* .11;
pdev = 0;
gap_int = 6;
pha_init = 0;
per_init = 600;

stim_seq = [repmat(stim_ioi,4,1)' stim_ioi*gap_int stim_comparison];

Wpha = 0.8;
Wper = 0.99;
Wpha_dec = Wpha;
Wper_dec = Wper;
WdecayConstant = .25;

pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;

chist = NaN(length(stim_seq)+1,1)';

M = table();
% model stim seq loop -------------------
s = 1; % stim counter
m = 1; % model counter
for i = 1:length(stim_seq)
    
    stim_vs_per = stim_seq(s)/per(m);
    nrep = floor(stim_vs_per);
    
    if nrep < 2 % normal model -- stim corresponds model
        
        C = mod(mod(pha(m),1) + mod(stim_vs_per,1),1);
        
        if C > .5
            C = C - 1;
        end
        
        % normal model equations
        pha(m+1) = (1 - Wpha)*C;
        per(m+1) = (1 + Wper*C)*per(m);
        chist(m) = C;
        corr_ok = 1;
        
        
        M.stim_index(m) = s;
        M.model_index(m) = m;
        M.stim_IOI(m) = stim_seq(s);
        M.stim_IOI_parts(m) = NaN;
        M.model_phase(m) = pha(m);
        M.model_period(m) = per(m);
        M.C(m) = C;
        M.was_corrected(m+1) = corr_ok;
        M.Wpha(m) = Wpha;
        M.Wper(m) = Wper;
        
        m = m+1;    % next model index
        
    else % if stim interval multiplies model interval (per(i))
        
        stim_int_parts = [repmat(per(m),1,nrep),  stim_seq(s)-per(m)*(nrep)];
        if stim_int_parts(end) == 0
            maxj = nrep;
        else
            maxj = nrep+1;
        end
        for j = 1:maxj % remainder from floor() func
            
            % chop long stim interval to exact multiples (and remainder)
            C = mod(mod(pha(m),1) + mod(stim_int_parts(j)/per(m),1),1);
            
            if C > .5
                C = C - 1;
            end
            

                % correction decay
                Wpha_dec = Wpha_dec - Wpha_dec * WdecayConstant;
                Wper_dec = Wper_dec - Wper_dec * WdecayConstant;
                M.Wpha(m) = Wpha_dec;
                M.Wper(m) = Wper_dec;
                
                pha(m+1) = (1 - Wpha_dec)*C;
                per(m+1) = (1 + Wper_dec *C)*per(m);
                chist(m) = C;
                corr_ok = -1;

            
            
            M.stim_index(m) = s;
            M.model_index(m) = m;
            if j == 1
                M.stim_IOI(m) = stim_seq(s);
            else
                M.stim_IOI(m) = NaN;
            end
            M.stim_IOI_parts(m) = stim_int_parts(j);
            M.model_phase(m) = pha(m);
            M.model_period(m) = per(m);
            M.C(m) = C;
            
            M.was_corrected(m+1) = corr_ok;
            
            m = m+1;
        end
        
    end
    
    s = s+1;
end

M = M(1:end-1,:);


corrected_times = logical([0 M.was_corrected(1:end)']);
stim_seq_times = [0 cumsum(stim_seq)];
per_seq_times = [0 cumsum(M.model_period)'];

f = figure;

subplot(2,1,1)
% oooooooooo STIM ooooooooooo
plot(stim_seq_times, ones(length(stim_seq_times),1) , '.k', 'MarkerSize', 30)
hold on;

% oooooooooo MODEL ooooooooooo
plot(per_seq_times, ones(length(per_seq_times),1)+.5 , '.r', 'MarkerSize', 30); hold on;
px = plot(per_seq_times(corrected_times), ones(length(per_seq_times(corrected_times)),1)+.5 , 'ok', 'MarkerSize', 20);
legend(px,'corrected times');

xticks(stim_seq_times); grid on;

yticks([1 1.5]); yticklabels({'stim', 'model'});
% set(gca, 'YGrid', 'off', 'xlim', [-100 max(stim_seq_times)+100]);
set(gca,'XTickLabelRotation',45)
EK_plotlabels('timepoints', '', ['C_{final} = ' num2str(C)], 20);
set(gcf,'Position',[1500 500 500 150])
xax = get(gca,'xlim');

subplot(2,1,2)
plot(per_seq_times,[NaN; NaN; M.Wpha(1:end-1)], '-o', 'linewidth', 3); hold on; 
plot(per_seq_times,[NaN; NaN; M.Wper(1:end-1)], '-o', 'linewidth', 3);
legend('Wpha (phase)', 'Wper (period)');
EK_plotlabels('timepoints', '', 'phase & period correction weights', 20);
set(gca, 'xlim', xax, 'ylim', [-0.1 1.1]);
xticks(stim_seq_times);

set(f,'Position',[1500         279         810         494])

%% bare model with modular stuff 

clear; clc; close all;


stim_ioi = 500;
stim_comparison = stim_ioi + stim_ioi* .11; % late 
pdev = 0;
gap_int = 1;
pha_init = 0;
per_init = 400;

stim_seq = [repmat(stim_ioi,4,1)' stim_ioi*gap_int stim_comparison];

Wpha = 0.8;
Wper = 0.99;
Wpha_dec = Wpha;
Wper_dec = Wper;

pha = NaN(length(stim_seq),1)'; pha(1) = pha_init;
per = NaN(length(stim_seq),1)'; per(1) = per_init;

for i = 1:length(stim_seq)-1
    
    C = mod((pha(i) + stim_seq(i)/per(i)),1);
    
    % [Q,R] = quorem(A,B) 
    
    if C > .5
        C = C - 1;
    end
    
    pha(i+1) = (1 - Wpha)*C;
    per(i+1) = (1 + Wper*C)*per(i);
    
end

% C_fin = C
enddiff = sum(per) - sum(stim_seq)
