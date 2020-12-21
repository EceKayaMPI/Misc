
clear,clc;
warning('off','all');
%% params

% Wphi = .5;              % 0 = no phase correction, 1 = full phase reset
% Wp = .4;                % 0 = no period corr. 1 = complete period corr where P matches curr IOI

all_models = table();
mdl = 1;

% Wphi_Wp = [0 1 .5 0 1 .5 0 1 .5; 0 0 0 1 1 1 .5 .5 .5]';

Wphi_Wp = [.5 .5 ; 0 1];

for par = 1:length(Wphi_Wp)
    
    Wphi = Wphi_Wp(par,1);
    Wp =  Wphi_Wp(par,2);
    
    all_models.Wphi(mdl:mdl+2) = Wphi;
    all_models.Wp(mdl:mdl+2) = Wp;
    
    
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
%                 condits.linc(c) = 1 - (.5 - abs(C));
            else
                condits.corr_resp(c) = 0;
            end
            
        elseif strcmp(condits.cmp_len{c}, 'comp same')
            if abs(C) < .01
                condits.corr_resp(c) = 1;
%                 condits.linc(c) = 1 - abs(C);
            else
                condits.corr_resp(c) = 0;
            end
            
        elseif strcmp(condits.cmp_len{c}, 'comp longer')
            if C > 0
                condits.corr_resp(c) = 1;
%                 condits.linc(c) = 1 - (.5 - abs(C));
            else
                condits.corr_resp(c) = 0;
            end
        end
        
    end
    apx = [];


%% linear trans

cfin = condits{:,7};

for c = 1:27
    if strcmp(condits.cmp_len{c}, 'comp shorter')
        condits.linc(c) = 1 - (.5 - abs(condits.C_fin(c)));
    elseif strcmp(condits.cmp_len{c}, 'comp same')
        condits.linc(c) = 1 - abs(condits.C_fin(c));
    elseif strcmp(condits.cmp_len{c}, 'comp longer')
        condits.linc(c) = 1 - (.5 - abs(condits.C_fin(c)));
    end
end

for c = 1:27

end



gamma = 4.5;

% a = softmax(n) = exp(n)/sum(exp(n))



% for i = 1:3:27
% 
%      denomsum = sum(exp(gamma*condits.linc(i:i+2)));
%      for k = 1:3
%      condits.luce(i+k-1) = (gamma*condits.linc(i+k-1))/denomsum;
%      end
% 
% end



tmp_beg = groupsummary(condits, {'cmp_beg', 'cmp_len'},'sum', {'linc'} );
lincgrp_beg = tmp_beg{:,3};

tmp_end = groupsummary(condits, {'st_end', 'cmp_len'},'sum', {'linc'} );
lincgrp_end = tmp_end{:,3};


for i = 1:3:height(tmp_beg)
    
    for k = 0:2
       tmp_beg.luce(i+k) = sum(tmp_beg.sum_linc(i:i+2));
       tmp_end.luce(i+k) = sum(tmp_end.sum_linc(i:i+2));
    end
end



% luce = table();
% for i = 1:3
%     luce.beginning(i) = lincgrp_beg(i)/sum(lincgrp_beg);
%     luce.ending(i) = lincgrp_end(i)/sum(lincgrp_end);
% end
% all_models.luce_beg(mdl:mdl+2) = luce.beginning;
% all_models.luce_end(mdl:mdl+2) = luce.ending;

%% plot

tittex = sprintf('Phase corr = %.2f, Period corr = %.2f', Wphi, Wp);

% figure;
% sgtitle(tittex);
% subplot(1,2,1)
% bar(luce.ending);
% xticklabels({'early','on-time','late'});
% xtickangle(45);
% set(gca, 'FontSize', 14);
% 
% subplot(1,2,2)
% bar(luce.beginning);
% xticklabels({'early','on-time','late'});
% xtickangle(45);
% set(gca, 'FontSize', 14);
% 
% 

% beg_manipulation = groupsummary(condits, {'cmp_beg'},'mean', {'linc'} );
% end_manipulation = groupsummary(condits, {'st_end'},'mean', {'linc'} );
% 
% figure;
% set(gca, 'FontSize',16);
% sgtitle(tittex);
% subplot(1,2,1)
% bar(end_manipulation.mean_linc);
% xticklabels(table2cell(end_manipulation(:,1)));
% xtickangle(45); %yline(10);
% set(gca, 'FontSize', 14);
% 
% subplot(1,2,2)
% bar(beg_manipulation.mean_linc);
% xticklabels(table2cell(beg_manipulation(:,1)));
% xtickangle(45); %yline(10);
% set(gca, 'FontSize', 14);

mdl = mdl+3;
end

%% plot all models together 

% legtext = cell(length(Wphi_Wp),1);
% 
% k = 1;
% for i = 1:3:height(all_models)
%     
%     legtext{k} = sprintf('Phase corr = %.2f, Period corr = %.2f', all_models{i,1}, all_models{i,2});
%     beggg(:,k) = all_models{i:i+2,3};
%     enddd(:,k) = all_models{i:i+2,4};
%     k = k+1;
% end
% 
% figure;
% set(gca, 'FontSize',16);
% subplot(1,2,1)
% bar(beggg);
% xticklabels({'early','on-time','late'});
% xtickangle(45);
% set(gca, 'FontSize', 14);
% title('beginning manipulation');
% 
% subplot(1,2,2)
% bar(enddd);
% xticklabels({'early','on-time','late'});
% xtickangle(45); 
% set(gca, 'FontSize', 14);
% title('ending manipulation');
% legend(legtext);

