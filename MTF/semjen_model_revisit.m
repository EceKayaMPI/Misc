clear, clc;

% m = params(1)         25
% t = params(2);        100
% alpha = params(3);    0.4
% beta = params(4);     0.3
real_params = [25, 100, 0.4, 0.3];

% Simulate data (100 experiments) with real params
all_ACVF = [];
for exp = 1:100 
    avg_exps = model(real_params,10);
    all_ACVF(exp,:) = calculate_ACVF(model(real_params,1)); % A_avg = model(params, nSeq)
end
% Compute real ACVF
real_ACVF = mean(all_ACVF,1);

%% try model with different params & find the one with the closest ACVF

% Define starting values and call parameter-estimation function
start_params = [30,100,1,1];
[final_params, final_discrepancy] = wrapper4fmin (start_params, real_ACVF);

%% functions

function [x, fval] = wrapper4fmin (start_params, real_ACVF)

[x, fval] = fminsearch(@compare_ACVF, start_params);

    function discrepancy = compare_ACVF(params)
        % model output (mean async. from nSeq sequences) with given params
        % (try diff. params here)
        model_out = model(params,10);
        % calculate ACVF of the model output with given params
        predicted_ACVF = calculate_ACVF(model_out);
        % difference between the real ACVF vs. predicted ACVF
        discrepancy = sum((predicted_ACVF - real_ACVF).^2);
    end
end

function A_avg = model(params, nSeq)
m = params(1);          % motor var (5^2)
t = params(2);          % timekeeper var (15^2)
alpha = params(3);      % 0.4
% beta = params(4);       % 0.3

N = 31;
k = 4; 
mDist = makedist('Gamma', 'a', k, 'b', sqrt(m/k));
tDist = makedist('Gamma', 'a', k, 'b', sqrt(t/k));

all_seqs = [];
for seq = 1:nSeq
    

    T = random(tDist, [N-1 1]);
    M = random(mDist, [N 1]);
    metronome = mean(T);
    
    for n = 1:30
        if n == 1
            % First Asynchrony is Metronome Interval
            A(n) = -metronome;
            % First-Order Error Correction
            Tx(n) = T(n) - alpha * A(n);
        else
            % Asynchrony of Previous Interval
            A(n) = sum(I(1:n-1)) - metronome*(n-1);
            % Second-Order Error Correction
            Tx(n) = T(n) - alpha * A(n);
        end
        % Intervals
        I(n) = Tx(n) + M(n+1) - M(n);
        
    end % n, taps
    I = I'; A = A';
    
    all_seqs(seq, :) = A;
    
end % experiment simulation

A_avg = mean(all_seqs,1);
end % sim function

function ACVF = calculate_ACVF(arr)
ACVF = [];
for lag = 0:3
    covmat = cov(arr(1:end-lag),arr(1+lag:end));
    ACVF(lag+1) = covmat(1,2);
end
end