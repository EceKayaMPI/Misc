%% Distribution Parameters
k = 4; m = 5^2; t = 15^2;
mDist = makedist('Gamma', 'a', k, 'b', sqrt(m/k));
tDist = makedist('Gamma', 'a', k, 'b', sqrt(t/k));

%% Experiment Parameters
N = 31; Ns = 200;
[I, A, Tx] = deal(zeros(N-1,1));

%% Synchronisation Parameters
metronome = mean(tDist);
alpha = .2:.05:1.8; beta = [-.2 -.1 0 .1 .2];
params = struct();

%% Big Chungus
idx = 1;
for a = 1:length(alpha)
    for b = 1:length(beta)
        [varA, corr1A, corr2A, varI, corr1I, corr2I] = deal(zeros(Ns,1));
        
        %% Simulate an Experiment
        for s = 1:Ns
            %% Sample from Distributions
            T = random(tDist, [N-1 1]);
            M = random(mDist, [N 1]);
            
            %% Generate Taps
            for n = 1:N-1
                % First Tap
                if n == 1
                    % First Asynchrony is Metronome Interval
                    A(n) = -metronome;
                    % First-Order Error Correction
                    Tx(n) = T(n) - alpha(a)*A(n);
                    
                % Subsequent Taps
                else
                    % Asynchrony of Previous Interval
                    A(n) = sum(I(1:n-1)) - metronome*(n-1);
                    % Second-Order Error Correction
                    Tx(n) = T(n) - alpha(a)*A(n) - beta(b)*A(n-1);
                end
                
                % Intervals
                I(n) = Tx(n) + M(n+1) - M(n);
            end
            
            varA(s) = var(A);
            corr1A(s) = corr(A(2:end), A(1:end-1));
            corr2A(s) = corr(A(3:end), A(1:end-2));

            varI(s) = var(I);
            corr1I(s) = corr(I(2:end), I(1:end-1));
            corr2I(s) = corr(I(3:end), I(1:end-2));
        end
        params(idx).alpha = alpha(a);
        params(idx).beta = beta(b);
        
        params(idx).varA = mean(varA);
        params(idx).corr1A = mean(corr1A);
        params(idx).corr2A = mean(corr2A);
        
        params(idx).varI = mean(varI);        
        params(idx).corr1I = mean(corr1I);
        params(idx).corr2I = mean(corr2I);
        
        idx = idx + 1;
    end
end

%% Graphs
estimator = fieldnames(params);
xlabs = {'var(A_{n})', 'p(A_{n}, A_{n+1})', 'p(A_{n}, A_{n+2})', 'var(I_{n})', 'p(I_{n}, I_{n+1})', 'p(I_{n}, I_{n+2})'};
for p = 1:6
    % Plot
    subplot(2,3,p);
    yline(0); hold on
    for i = 1:5
        plot([params(i:5:end).alpha], [params(i:5:end).(estimator{p+2})]); 
    end
    
    % Axes
    xlim([.2 1.8]);
    if mod(p-1, 3) == 0
        ylim([100 1450]);
    else
        ylim([-1 1]);
    end
    xlabel(xlabs{p}, 'Interpreter', 'tex');
end