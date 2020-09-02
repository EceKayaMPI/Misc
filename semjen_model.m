
% m = params(1)         25
% t = params(2);        100
% alpha = params(3);    0.4    
% beta = params(4);     0.3
real_params = [25, 100, 0.4, 0.3];

% Simulate data & Compute real ACVF 
real_ACVF = model (real_params);

% Define starting values and call parameter-estimation function
start_params = [1,1,1,1];
[final_params, final_discrepancy] = wrapper4fmin (start_params, real_ACVF);

function [x, fval] = wrapper4fmin (start_params, real_ACVF)

[x, fval] = fminsearch(@bof, start_params);

    function rmsd = bof(params)
        predicted_ACVF = model (params);
        sd = (predicted_ACVF - real_ACVF).^2;
        rmsd = sqrt(sum(sd)/numel(sd));
    end
end

function ACVF = model (params)
m = params(1);
t = params(2);
alpha = params(3);
beta = params(4);

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