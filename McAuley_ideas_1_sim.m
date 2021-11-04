
cd /Users/ece.kaya/MATLAB/Githubbed/Misc

stim_ioi = 500;
pdev = .11;

gap_int = 6;
pha_init = 0;
per_init = 600;

stim_comparison = stim_ioi + stim_ioi* pdev;
stim_seq = [repmat(stim_ioi,4,1)' stim_ioi*gap_int stim_comparison];

Wpha = 0.8;
Wper = 0.99;

C_fin = McAuley_freezer (stim_seq, pha_init, per_init, Wpha, Wper);

%% simulate exp & model resp
load /Users/ece.kaya/MATLAB/Githubbed/OscPerc/OPC_matlab/OPC_50.mat

ioilist = OPC_50(1).S(2).condits(:,1).*1000;
pdevlist = OPC_50(1).S(2).condits(:,2);
condits = table(ioilist,pdevlist);

condits = condits(logical(ioilist<700) & logical(ioilist>300),:);

sim = table();

% fixed params
gap_int = 6;
pha_init = 0;
per_init = 500;

% fixed for simulation
Wpha = 0.7;
Wper = 0.99;

for i = 1:height(condits)
    
    stim_ioi = condits.ioilist(i);
    pdev = condits.pdevlist(i);
    
    stim_comparison = stim_ioi + stim_ioi* pdev;
    stim_seq = [repmat(stim_ioi,4,1)' stim_ioi*gap_int stim_comparison];
    
    C_fin = McAuley_freezer (stim_seq, pha_init, per_init, Wpha, Wper);
    
    sim.ioi(i) = stim_ioi;
    sim.dev(i) = pdev;
    
    
    sim.cfin(i) = C_fin;
    sim.discrep(i) = C_fin - pdev;
    
    if sim.cfin(i) > 0      % if Cfinal is positive 
        sim.resp(i) = 1;   % model will think the interval is 'shorter' than predicted (?)
    elseif sim.cfin(i) < 0  % if Cfinal negative 
        sim.resp(i) = -1;    
    else
        sim.resp(i) = NaN;
    end
    
    
    
    % resp coding
    if pdev > 0
        sim.devdir(i) = 1;
        if sim.resp(i) > 0
            sim.hit(i) = 1;
            sim.bias(i) = 0;
        elseif sim.resp(i) < 0
            sim.hit(i) = 0;
            sim.bias(i) = -1;
        end
    elseif pdev < 0
        sim.devdir(i) = -1;
         if sim.resp(i) > 0
            sim.hit(i) = 0;
            sim.bias(i) = 1;
        elseif sim.resp(i) < 0
            sim.hit(i) = 1;
            sim.bias(i) = 0;
        end       
        
    end
    
end    

% ---------------- HIT MEAS ------------------------
ws = 30;
[hitmm, ioimm] = EK_movmean_by(sim.hit,sim.ioi,ws);
[hitmmpos, ioimmpos] = EK_movmean_by(sim.hit(sim.devdir==1), sim.ioi(sim.devdir==1),ws);
[hitmmneg, ioimmneg] = EK_movmean_by(sim.hit(sim.devdir==-1), sim.ioi(sim.devdir==-1),ws);

figure;
sim = sortrows(sim,'ioi','ascend');
plot(sim.ioi, sim.discrep);


figure;
sgtitle({'Model simulation responses', ['initial period = ' num2str(per_init)],...
    ['phase correction = ' num2str(Wpha)],  ['period correction = ' num2str(Wper)]}, 'FontSize', 22);
subplot(3,1,1)

plot(ioimm, hitmm, 'linewidth', 2); hold on;
yline(0.5, '--k', 'linewidth', 2);
EK_plotlabels('IOI', 'hit movmean (ws=30)', 'all resp ', 16);

subplot(3,1,2)
plot(ioimmpos, hitmmpos, 'linewidth', 2); hold on;
plot(ioimmneg, hitmmneg, 'linewidth', 2);
legend('positive dev', 'negative dev');
EK_plotlabels('IOI', 'hit movmean (ws=30)', 'pos-neg resp ', 16);

% ---------------- BIAS MEAS ------------------------
[biasmm, ioimm] = EK_movmean_by(sim.bias, sim.ioi,ws);
absbiasmm = abs(biasmm);

subplot(3,1,3)
plot(ioimm, biasmm, 'linewidth', 2); hold on;
yline(0, '--k', 'linewidth', 2);
EK_plotlabels('IOI', 'bias movmean (ws=30)', 'all resp ', 16);
    
    
    %% the function
    
    function C_fin = McAuley_freezer (stim_seq, pha_init, per_init, Wpha, Wper)
    pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
    per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;
    
    chist = NaN(length(stim_seq)+1,1)';
    
    % M = table();
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
            
            
            %         M.stim_index(m) = s;
            %         M.model_index(m) = m;
            %         M.stim_IOI(m) = stim_seq(s);
            %         M.stim_IOI_parts(m) = NaN;
            %         M.model_phase(m) = pha(m);
            %         M.model_period(m) = per(m);
            %         M.C(m) = C;
            %         M.was_corrected(m+1) = corr_ok;
            
            
            
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
                
                if j < maxj
                    % model doesn't correct for phase & period here (Wper = 0, Wpha =
                    % 0)
                    pha(m+1) = (1 - 0)*C;
                    per(m+1) = (1 + 0 *C)*per(m);
                    chist(m) = C;
                    corr_ok = 0;
                else
                    % normal model equations to predict final interval (bcs
                    % there will be stim next )
                    pha(m+1) = (1 - Wpha)*C;
                    per(m+1) = (1 + Wper*C)*per(m);
                    chist(m) = C;
                    corr_ok = 1;
                    
                end
                
                
                %             M.stim_index(m) = s;
                %             M.model_index(m) = m;
                %             if j == 1
                %                 M.stim_IOI(m) = stim_seq(s);
                %             else
                %                 M.stim_IOI(m) = NaN;
                %             end
                %             M.stim_IOI_parts(m) = stim_int_parts(j);
                %             M.model_phase(m) = pha(m);
                %             M.model_period(m) = per(m);
                %             M.C(m) = C;
                %             M.was_corrected(m+1) = corr_ok;
                
                m = m+1;
            end
            
        end
        
        s = s+1;
    end
    
    C_fin = C;
    
    end
    
    
    function C_fin = McAuley_main (stim_seq, pha_init, per_init, Wpha, Wper)
    pha = NaN(length(stim_seq)+1,1)'; pha(1) = pha_init;
    per = NaN(length(stim_seq)+1,1)'; per(1) = per_init;
    
    chist = NaN(length(stim_seq)+1,1)';
    
    for i = 1:length(stim_seq)

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
            
            
            %         M.stim_index(m) = s;
            %         M.model_index(m) = m;
            %         M.stim_IOI(m) = stim_seq(s);
            %         M.stim_IOI_parts(m) = NaN;
            %         M.model_phase(m) = pha(m);
            %         M.model_period(m) = per(m);
            %         M.C(m) = C;
            %         M.was_corrected(m+1) = corr_ok;
            
            
            
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
                
                if j < maxj
                    % model doesn't correct for phase & period here (Wper = 0, Wpha =
                    % 0)
                    pha(m+1) = (1 - 0)*C;
                    per(m+1) = (1 + 0 *C)*per(m);
                    chist(m) = C;
                    corr_ok = 0;
                else
                    % normal model equations to predict final interval (bcs
                    % there will be stim next )
                    pha(m+1) = (1 - Wpha)*C;
                    per(m+1) = (1 + Wper*C)*per(m);
                    chist(m) = C;
                    corr_ok = 1;
                    
                end
                
                
                %             M.stim_index(m) = s;
                %             M.model_index(m) = m;
                %             if j == 1
                %                 M.stim_IOI(m) = stim_seq(s);
                %             else
                %                 M.stim_IOI(m) = NaN;
                %             end
                %             M.stim_IOI_parts(m) = stim_int_parts(j);
                %             M.model_phase(m) = pha(m);
                %             M.model_period(m) = per(m);
                %             M.C(m) = C;
                %             M.was_corrected(m+1) = corr_ok;
                
                m = m+1;
            end
            
        end
        
        s = s+1;
    end
    
    C_fin = C;
    
    end