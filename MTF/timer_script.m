close all; clear; clc;

kappa = .1;

n = 3;      % number of osc
times = [];
phases = [];
waves = [];

for i = 1:n
    phases(1,i) = 2*pi*rand(1);
end


figure;
set(gcf,'Position',[100 100 500 1000])
subplot(2,1,1)
k = 1;
for t = 1:.2:10 % time
    
    for i = 1:n
        if k == 1
            pha = phases(k,i);
            
        else
            
            sum_sin_pha_diff = 0;
            for j = [1:i-1 i+1:n]
                sum_sin_pha_diff = sum_sin_pha_diff + ...
                    sin(phases(k-1,j) - phases(k-1,i));
            end
            pha = phases(k-1,i) + kappa * sum_sin_pha_diff;
            
        end
        
        times(k) = t;
        phases(k,i) = pha;
        waves(k,i) = sin(t + pha);
        
        
        make_dots(t + pha); 
        drawnow
        hold on;
        
    end
    
    
    hold off;
    k = k+1;
end

subplot(2,1,2)
for i = 1:n
    plot(times, waves(:,i)); hold on;
end




%%
close all; clear; clc;

kappa = .01;
osc = table();

figure;
subplot(2,2,1)
k = 1;
for t = 1:.5:50 % sampling rate
    
    osc.t(k) = t;
    
    if k == 1
        pha_1 = 0;
        pha_2 = 1;
        
    else
        pha_1 = osc.pha_1(k-1) + kappa * sin(osc.pha_2(k-1) - osc.pha_1(k-1));
        pha_2 = osc.pha_2(k-1) + kappa * sin(osc.pha_1(k-1) - osc.pha_2(k-1));
        
    end
    
    osc.pha_1(k) = pha_1;
    osc.pha_2(k) = pha_2;
    
    osc.wave_1(k) = sin(t + pha_1);
    osc.wave_2(k) = sin(t + pha_2);
    
    make_dots(pha_1); hold on;
    make_dots(pha_2); hold off;
    drawnow
    
    k = k+1;
end

subplot(2,2,2)
plot(osc.t, osc.wave_1); hold on;
plot(osc.t, osc.wave_2);


%%
close all; clear; clc;

kappa = .5;
osc = table();

figure;
subplot(2,2,1)
k = 1;
for t = 1:.1:10 % sampling rate
    
    pha_1 = t;
    pha_2 = t+1;
    
    osc.t(k) = t;
    osc.theta_1(k) = pha_1;
    osc.theta_2(k) = pha_2;
    
    
    
    
    make_dots(pha_1); hold on;
    make_dots(pha_2); hold off;
    drawnow
    
    k = k+1;
end

subplot(2,2,2)
plot(osc.t, sin(osc.theta_1)); hold on;
plot(osc.t, sin(osc.theta_2));




%% functions

% function kuramoto(thetaVec, k,

function make_dots(theta)
plot(cos(theta), sin(theta), "Marker",".","MarkerSize",50);
set(gca,'xlim', [-1 1], 'ylim', [-1 1]);
end


