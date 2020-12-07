
clear,clc;
%% params

Wphi = .5;              % 0 = no phase correction, 1 = full phase reset
Wp = .4;                % 0 = no period corr. 1 = complete period corr where P matches curr IOI

pha = 0;            % phase corr term
per = 600;          % period corr term
C = 0;              % temporal contrast

%% appendix B table

condits = table();

contextIOI = 600;
context_seq = repmat(contextIOI,1,7);

standard_opts = [540, 600, 660];

k  = 1;

for i = 1:3 % 3 comparison beg (early/on-time/late)
    
    for s = 1:3 % 3 standard ending
        standard = standard_opts(s);            % select standard
        ISI = 1200 - (standard - contextIOI);   % compensate ISI for this standard
        
        
        for c = 1:3 % 3 x comparison length     % manipulate comp length
            
            if c == 1 % comp shorter
                condits.cmp_len{k} = 'comp shorter';
                comparison = standard - standard*.1;
                
            elseif c == 2 % comp same
                
                condits.cmp_len{k} = 'comp same';
                comparison = standard;
                
            else % comp longer
                condits.cmp_len{k} = 'comp longer';
                comparison = standard + standard*.1;
            end
            
            
            if i == 1
                condits.cmp_beg{k} = 'early beginning';
                ISI_fin = ISI - abs(comparison - standard);
                
            elseif i == 2
                condits.cmp_beg{k} = 'on-time beginning';
                ISI_fin = ISI;
                
            elseif i == 3
                condits.cmp_beg{k} = 'late beginning';
                ISI_fin = ISI + abs(comparison - standard);
            end
            
            if s == 1
                condits.st_end{k} = 'early ending';
                
            elseif s == 2
                condits.st_end{k} = 'on-time ending';
                
            elseif s == 3
                condits.st_end{k} = 'late ending';
            end
            
            condits.st(k) = standard;
            condits.ISI(k) = ISI_fin;
            condits.com(k) = comparison;
            
            
            grp = [standard, ISI_fin, comparison];
            apx(:,k) = grp;
            k = k+1;
        end
    end
end
apx = apx.';

%% model & find Cfin

for c = 1:27
    
    stim_seq = [context_seq apx(c,:)];
    
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
    
    condits.C_fin(c) = C;
    
    
    if strcmp(condits.cmp_len{c}, 'comp shorter')
        if C < 0
            condits.corr_resp(c) = 1;
        else
            condits.corr_resp(c) = 0;
        end
    elseif strcmp(condits.cmp_len{c}, 'comp same')
        if C > 0
            condits.corr_resp(c) = 1;
        else
            condits.corr_resp(c) = 0;
        end
    end
    
end

beg_manipulation = groupsummary(condits, {'cmp_beg'},'sum', {'corr_resp'} );
end_manipulation = groupsummary(condits, {'st_end'},'sum', {'corr_resp'} );

% plot
figure;
subplot(1,2,1)
bar(end_manipulation.sum_corr_resp); 
xticklabels(table2cell(end_manipulation(:,1)));
yline(10);set(gca, 'FontSize', 14);

subplot(1,2,2)
bar(beg_manipulation.sum_corr_resp); 
xticklabels(table2cell(beg_manipulation(:,1)));
yline(10);set(gca, 'FontSize', 14);

