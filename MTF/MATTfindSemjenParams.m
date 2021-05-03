paramsToFind = [.4 25 100];
% startValues = [2 1 25 49];
startValues = [.3 20 100];
acvf = SemjenFunction(paramsToFind, 2000);

estimate = paramSearch(acvf, startValues);

function out = SemjenFunction(in, Ns)
%% Distribution Parameters
k = 4;
mDist = makedist('Gamma', 'a', k, 'b', sqrt(in(2)/k));
tDist = makedist('Gamma', 'a', k, 'b', sqrt(in(3)/k));

%% Experiment Parameters
N = 31; acvf = zeros(Ns,4);
metronome = mean(tDist);
[I, A, Tx] = deal(zeros(N-1,1));

%% Simulate Experiments
for s = 1:Ns
    % Sample from Distributions
    T = random(tDist, [N-1 1]);
    M = random(mDist, [N 1]);

    % Generate Taps
    for n = 1:N-1
        % First Tap
        if n == 1
            % First Asynchrony is Metronome Interval
            A(n) = -metronome;
            % First-Order Error Correction
            Tx(n) = T(n) - in(1)*A(n);
            
        % Subsequent Taps
        else
            % Asynchrony of Previous Interval
            A(n) = sum(I(1:n-1)) - metronome*(n-1);
            % Second-Order Error Correction
            Tx(n) = T(n) - in(1)*A(n); % - in(2)*A(n-1);
        end
        
        % Intervals
        I(n) = Tx(n) + M(n+1) - M(n);
    end
    acvf(s,:) = ACVF(A);
end
out = mean(acvf);
end

%% func
function out = paramSearch(acvf, x0)
f = @(params,acvf)sum((SemjenFunction(params,100)-acvf).^2);
fun = @(x)f(x,acvf);

options = optimset(...
    'Display', 'iter',...
    'PlotFcns', @optimplotfval,...
    'MaxFunEvals', 3000,...
    'MaxIter', 1000);
out = fminsearch(fun, x0, options);

% out = fminsearch(fun, x0)

end

function ACVF = ACVF(arr)
ACVF = [];
for lag = 0:3
    covmat = cov(arr(1:end-lag),arr(1+lag:end));
    ACVF(lag+1) = covmat(1,2);
end
end