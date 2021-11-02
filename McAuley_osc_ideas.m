

%% oscillator free during intervals that double its period
clear; clc; close all;

fs = 44100;
Ts = 1/fs;
sinwa = [];

stim_ioi = 500;
gap_int = 6;
pha_init = 0;
per_init = 600;

stim_seq = [repmat(stim_ioi,4,1)' stim_ioi*gap_int stim_ioi];

Wpha = 0.4;
Wper = 0.99;

pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;

chist = NaN(length(stim_seq)+1,1)';


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
        
        %         s = s+1;    % next stim index
        m = m+1;    % next model index
        
    else % if stim interval multiplies model interval (per(i))
        
        stim_int_parts = [repmat(per(m),1,nrep-1),  stim_seq(s)-per(m)*(nrep-1)];
        for j = 1:nrep % model oscillates alone that many times
            
            % chop long stim interval to exact multiples (and remainder)
            C = mod(mod(pha(m),1) + mod(stim_int_parts(j)/per(m),1),1);
            
         if C > .5
            C = C - 1;
        end           
            
            
            % model doesn't correct for phase & period here (Wper = 0, Wpha =
            % 0)
            pha(m+1) = (1 - 0)*C;
            per(m+1) = (1 + 0 *C)*per(m);
            chist(m) = C;
            
            m = m+1;
        end
        
    end
    
    s = s+1;
end

stim_seq_times = [0 cumsum(stim_seq)];
per_seq_times = [0 cumsum(per)];

figure;
plot(stim_seq_times, ones(length(stim_seq_times),1) , '.k', 'MarkerSize', 30)
hold on;
% plot(sinwa, 'linewidth', 1.5)

plot(per_seq_times, ones(length(per_seq_times),1)+.5 , '.r', 'MarkerSize', 30)
set(gcf,'Position',[1000 500 900 50])

%% oscillator doubles its period when stim interval is double

clear; clc; close all;

fs = 44100;
Ts = 1/fs;
sinwa = [];

stim_ioi = 500;
gap_int = 5;
pha_init = 0;
per_init = 500;

stim_seq = [repmat(stim_ioi,5,1)' stim_ioi*gap_int stim_ioi];

Wpha = 0.56;
Wper = 0.52;

pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;

chist = NaN(length(stim_seq)+1,1)';

% model stim seq loop -------------------
for i = 1:6%length(stim_seq)
    
    if per(i) < stim_seq(i)
        nrep = floor(stim_seq(i)/per(i));

        if nrep > 1
            per(i) = per(i)*nrep;
        end
    elseif per(i) > stim_seq(i)
        nrep = floor(per(i)/stim_seq(i));

        if nrep > 1
            per(i) = per(i)/nrep;
        end
    end
    
    C = mod(mod(pha(i),1) + mod(stim_seq(i)/per(i),1),1);
    
    if C > .5
        C = C - 1;
    end
    
    % normal model equations
    pha(i+1) = (1 - Wpha)*C;
    per(i+1) = (1 + Wper*C)*per(i); % << separate currnt period to use vs next to calc
    chist(i) = C;
    
    
    
    
    
    
    
    
end

stim_seq_times = [0 cumsum(stim_seq)];
per_seq_times = [0 cumsum(per)];

figure;
plot(stim_seq_times, ones(length(stim_seq_times),1) , '.k', 'MarkerSize', 30)
hold on;
% plot(sinwa, 'linewidth', 1.5)

plot(per_seq_times, ones(length(per_seq_times),1)+.5 , '.r', 'MarkerSize', 30)
set(gcf,'Position',[1000 500 900 50])

%% bare model

clear; clc; close all;

stim_ioi = 500;
gap_int = 5;
pha_init = 0;
per_init = 600;

stim_seq = [repmat(stim_ioi,5,1)' stim_ioi*gap_int stim_ioi];

Wpha = 0.99;
Wper = 0.99;

pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;

chist = NaN(length(stim_seq)+1,1)'; 

% model stim seq loop -------------------
for i = 1:length(stim_seq)

    
        C = mod((pha(i) + stim_seq(i)/per(i)),1);
    
        if C > .5
            C = C - 1;
        end
        
        % normal model equations
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
% set(gca, 'YGrid', 'off', 'xlim', [-100 max(stim_seq_times)+100]);
set(gca,'XTickLabelRotation',45)
EK_plotlabels('timepoints', '', ['C_{final} = ' num2str(C)], 15);
set(gcf,'Position',[1500 500 500 150])