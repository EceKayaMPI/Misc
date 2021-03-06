
%% simulate data with established params and get the ACVF to be matched

[simulated_I] = simulate_sync (25,100,0.4,0.3);
[ACVF_goal] = calculate_ACVF(simulated_I);


% [pval,fval] = fminsearch(@(params) sumsq_ACVF_diff(params(1),params(2),params(3),params(4)),[25, 100, 0.3, 0.3]);
[acvf_diff] = sumsq_ACVF_diff (25,100,0.4,0.3,ACVF_goal );
%% function to minimize

function [estim_params, discrepancy] = fminsearcher (ACVF_goal)

function [acvf_diff] = sumsq_ACVF_diff (m,t,alpha,beta)
[test_I] = simulate_sync (m,t,alpha,beta);
[test_ACVF] = calculate_ACVF(test_I);
% squared deviations between the two acvf arrays (summed sq deviations)
difference = ACVF_goal - test_ACVF;
acvf_diff = sum(difference(:).^2);
end

end
%% the function that simulates data only (out = I, in = m,t,alpha,beta)
function [I] = simulate_sync (m,t,alpha,beta)
% m = 25;
% t = 100;
% alpha = 0.4;
% beta = 0.3;

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
end
I = I'; A = A';
end
%% func calculates ACVF lags 0-3 (out = ACVF, in = I)
function [ACVF] = calculate_ACVF(I)
ACVF = [];
lag = 0;
while length(ACVF) < 4
    x0 = I(1:end-lag);
    x1 = I(1+lag:end);
    covmat = cov(x0,x1);
    ACVF(lag+1) = covmat(1,2);
    lag = lag +1;
end
end