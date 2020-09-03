clear, clc;

% m = params(1)         25
% t = params(2);        100
% alpha = params(3);    0.4    
% beta = params(4);     0.3
real_params = [25, 100, 0.4, 0.3];

% Simulate data with real params 
% Compute real ACVF 
all_ACVF = [];
for exp = 1:100
all_ACVF(exp,:) = calculate_ACVF(model(real_params,1)); % I_avg = model(params, nSeq)
end
real_ACVF = mean(all_ACVF,1);

% Define starting values and call parameter-estimation function
start_params = [50,50,1,1];
[final_params, final_discrepancy] = wrapper4fmin (start_params, real_ACVF);

function [x, fval] = wrapper4fmin (start_params, real_ACVF)

% options = optimset('Display','iter','PlotFcns',@optimplotfval);
[x, fval] = fminsearch(@real_vs_estimated, start_params);

    function discrepancy = real_vs_estimated(params)
        % calculate ACVF of the model output with given params
        predicted_ACVF = calculate_ACVF(model(params,1));
        % difference between the real ACVF vs. predicted ACVF
        discrepancy = sum((predicted_ACVF - real_ACVF).^2);
    end
end

function A_avg = model(params, nSeq)
m = params(1);
t = params(2);
alpha = params(3);
beta = params(4);

all_seqs = [];
for experiment = 1:nSeq
NI = 31;
k = 4;
M = gamrnd(k,sqrt(m/k),[NI 1]);
T = gamrnd(k,sqrt(t/k),[NI-1 1]);
metronome = mean(T);

    for n = 1:30
        if n == 1
            A(n) = gamrnd(k,sqrt(m/k),[1 1]); % some motor delay
            Tx(n) = T(n) -  alpha * A(n);
            I(n) = Tx(n) + M(n+1) - M(n);
        else
            A(n) = sum(I(1:n-1)) - metronome*(n-1);
            Tx(n) = T(n)  - alpha * A(n) - beta * A(n-1);
            I(n) = Tx(n) + M(n+1) - M(n);
        end
    end % n, taps
    I = I'; A = A';

all_seqs(experiment, :) = A;
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