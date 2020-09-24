
clear,clc;

N = 34;
stim_IOIs = linspace(600,100,N);

alpha =0.4;
beta =0.1;
m =0.5;
gamma =0.4;

stim = [0 cumsum(stim_IOIs)];

t_pred = zeros(1,N);

for n = 3:N

%         fit = polyfit(n-2:n,diff(stim(n-3:n)),1);     
%         IOI_pred(n+1) = fit(2) + fit(1) * (n+1);
%         IOI_track(n+1) = t(n) - t(n-1);
%         t_pred(n+1) = stim_times(n)+ IOI_pred(n+1);
%         t_pred(n+1) = t(n) + (m * IOI_pred(n+1) + (1-m) * IOI_track(n+1)) + TK2(n);
        
% basic version
        fit = polyfit(n-2:n,stim_IOIs(n-2:n),1);     
        IOI_pred(n+1) = fit(2) + fit(1) * (n+1)
        t_pred(n+1) = stim(n)+ IOI_pred(n+1)

end 

figure;
subplot(2,1,1)
plot(diff(stim)), hold on; plot(diff(t_pred)), hold off;
subplot(2,1,2)
EK_dotplot(stim); hold on; EK_dotplot(t_pred); legend('stimuli', 'ADAM predicted'); hold off;