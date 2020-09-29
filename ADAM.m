clear,clc;

N = 34;
stim_IOIs = linspace(1600,100,N);
stim = [0 cumsum(stim_IOIs)];

limits = [250, 1000];

alpha =0.4;
beta =0.1;
m =0.5;                 % 0 = track | 1 = predict 
gamma_init =0.2;        % 0 = anticipate | 1 = adapt 

mvarProp = .5;          % 0 = motor flexible | 1 = motor unflexible 

T = zeros(1,N);         % timekeeper (interval)
t = zeros(1,N);         % ADAM's tap (timestamp)
asyn = zeros(1,N);      % asynchrony (stimtime - t)    
t_pred = zeros(1,N);    % output tap time of anticipation module
t_fin = zeros(1,N);     % final output tap time from the joint module
IOI_pred = zeros(1,N);  % predicted IOI (via anticipation)
IOI_track = zeros(1,N); % tracked IOI (previous taps)

TK1 = normrnd(0,1,[1,N]);            % sampled from a normal dist
TK2 = normrnd(0,1,[1,N]);            % sampled from a normal dist
M = gamrnd(4,sqrt(25/4),[1,N+1]);      % drawn from gamma dist 
M_init = M;

T(2:3) = stim_IOIs(1:2) + TK1(1:2);  % initialize timekeeper with first 2 stims


for n = 2:N
    gamma = gamma_init;

    % adaptation ----------------------------------------
    asyn(n) = stim(n-1) - t(n-1);                     % most recent asynchrony
    
    t(n+1) = t(n) + ( T(n) + (alpha + beta) * asyn(n) + TK1(n));    % eq. (11)
    T(n+1) = T(n) + beta * asyn(n);                                 % eq. (12)

    if n > 2
    % anticipation --------------------------------------    
        fit = polyfit(n-2:n,stim_IOIs(n-2:n),1);     
        IOI_pred(n+1) = fit(2) + fit(1) * (n+1);
        IOI_track(n+1) = t(n) - t(n-1);

        t_pred(n+1) = t(n) + (m * IOI_pred(n+1) + (1-m) * IOI_track(n+1)) + TK2(n);
        
        % motor variance increases proportionally after predicted (next) IOI < min(limits)
        if IOI_pred(n+1) < min(limits)
            
            M(n+1) =  M(n+1) + mvarProp * (min(limits) - IOI_pred(n+1));

        end

        
    else % n = 2 not enough stim to predict future
        gamma = 1; 
    end
   
    % joint -------------------------------------- 
    ASYN(n+1) = t(n+1) - t_pred(n+1);
    t_fin(n+1) = t(n+1) - ((1-gamma) * ASYN(n+1)) + M(n) - M(n-1);

end

figure;
subplot(2,1,1)
plot(diff(stim)), hold on; plot(diff(t_fin)), hold off;
subplot(2,1,2)
EK_dotplot(stim); hold on; EK_dotplot(t_fin); legend('stimuli', 'ADAM'); hold off;