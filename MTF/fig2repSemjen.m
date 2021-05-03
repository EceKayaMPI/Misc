clear, clc;
% Histograms and scatter plots of parameter estimates in N = 1000
% simulated experiments (t = 100 m = 25 a = .4 b = .3)

all_params = zeros (1000,4);
for experiment = 1:1000

real_params = [25, 100, 0.4, 0.3];
% m = params(1)         25
% t = params(2);        100
% alpha = params(3);    0.4    
% beta = params(4);     0.3

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
% chose params
% final_params ~ params << should be. 


all_params(experiment,:) = final_params;
end




function [x, fval] = wrapper4fmin (start_params, real_ACVF)

% options = optimset('Display','iter','PlotFcns',@optimplotfval);
[x, fval] = fminsearch(@real_vs_estimated, start_params);

    function discrepancy = real_vs_estimated(params)
        % calculate ACVF of the model output with given params
        sim_ACVF = calculate_ACVF(model(params,1));
        % difference between the real ACVF vs. predicted ACVF
        discrepancy = sum((sim_ACVF - real_ACVF).^2);
    end
end

function I_avg = model(params, nSeq)
m = params(1);
t = params(2);
alpha = params(3);
beta = params(4);

all_seqs = [];
for seq = 1:nSeq
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

all_seqs(seq, :) = I;
end % experiment simulation
I_avg = mean(all_seqs,1);
end % sim function 

function ACVF = calculate_ACVF(I)
ACVF = [];
for lag = 0:3
    covmat = cov(I(1:end-lag),I(1+lag:end));
    ACVF(lag+1) = covmat(1,2);
end
end

