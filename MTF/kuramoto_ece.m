
%% stim??

close all; clear; clc;

kappa = .1;

maxtime = 30;
timeInc = .2;

nsamps = length(1:timeInc:maxtime);

% generate stim 
stim = [];
while 1
    
    stim = [stim ones(10,1)' zeros(20,1)'];
    
    if length(stim) > nsamps
        stim = stim(1:nsamps);
        break
    end
end


osc = table();
osc.pha_1 = stim';

figure;
subplot(2,1,1)
k = 1;
for t = 1:timeInc:maxtime % sampling rate
    
    osc.t(k) = t;
    
    if k == 1
        pha_1 = osc.pha_1(k);
        pha_2 = 0;
        
    else
        pha_1 = osc.pha_1(k);
        pha_2 = osc.pha_2(k-1) + kappa * sin(osc.pha_1(k-1) - osc.pha_2(k-1));
        
    end
    
    osc.pha_1(k) = pha_1;
    osc.pha_2(k) = pha_2;
    
    osc.wave_1(k) = pha_1;
    osc.wave_2(k) = sin(t + pha_2);
    
    make_dots(pha_1); hold on;
    make_dots(pha_2); hold off;
    drawnow
    
    k = k+1;
end

subplot(2,1,2)
plot(osc.t, osc.pha_1); hold on;
plot(osc.t, osc.wave_2);



%% kappas ???
close all; clear; clc;

clrs = lines(10);


K = .5;

n = 5;      % number of osc
times = [];
phases = [];
waves = [];

kappas = [0    0    0    0    0 ;...      % 1st osc (stim) is independent
          K    0    K    K    K;...      % 2nd (auditory) 
          K    K    0    K    K;...
          K    K    K    0    K;...          
          K    K    K    K    0]; 

      
% kappas = [0    0    0;...      % 1st osc (stim) is independent
%           K    0    K;...      % 2nd (auditory) 
%           K    K    0];       
      
      
% for i = 1:n
%     phases(1,i) = 2*pi*rand(1);
% end
phases = [(1:n)'/n*2*pi]';

figure;
set(gcf,'Position',[100 100 500 1000])
subplot(2,1,1)
k = 1;
for t = 1:.5:50 % time
    
    for i = 1:n
        if k == 1
            pha = phases(k,i);
            
        else
            
            sum_sin_pha_diff = 0;
            for j = 1:n
                sum_sin_pha_diff = sum_sin_pha_diff + ...
                kappas(j,i) * sin(phases(k-1,j) - phases(k-1,i));
            % ??? how to put kappa vector ?? 
            end
            
            pha = phases(k-1,i) + sum_sin_pha_diff;
            
        end
        
        times(k) = t;
        phases(k,i) = pha;
        waves(k,i) = sin(t + pha);
        
        
        make_dots(t + pha, clrs(i,:)); 
        drawnow
        hold on;
        
    end
    
    
    hold off;
    k = k+1;
end

subplot(2,1,2)
for i = 1:n
    plot(times, waves(:,i), 'Color', clrs(i,:)); hold on;
end






%% basic

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
function make_dots(theta, clr)
if nargin<2
    clr = [0 0 0];
end
plot(cos(theta), sin(theta), "Marker",".","MarkerSize",50, "MarkerFaceColor", clr);
set(gca,'xlim', [-1 1], 'ylim', [-1 1]);
end


