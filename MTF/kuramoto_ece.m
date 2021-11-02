%% stim weird

close all; clear; clc;
cd '/Users/ece.kaya/MATLAB/Githubbed/Misc/MTF';

kappa = .1;

maxtime = 33;
timeInc = .1;

nsamps = length(1:timeInc:maxtime);

% generate stim 
stim = [];
stimSam = 10;
silSam = 20;


while 1
    
    stim = [stim zeros(silSam,1)' ones(stimSam,1)' ];
    
    if sum(islocalmax(stim)) > 4
        stim = [stim zeros(silSam,1)' zeros(silSam,1)' zeros(silSam,1)' ones(stimSam,1)' zeros(silSam,1)'];
    end
    
    if length(stim) > nsamps
        stim = stim(1:nsamps);
        break
    end
end

figure
stem(stim);

osc = table();
osc.pha_1 = stim';

figure;
set(gcf,'Position',[100 100 500 1000]);


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
plot(osc.t, osc.pha_1, 'linewidth', 2); hold on;
plot(osc.t, osc.wave_2, 'linewidth', 2);
%% stim??

close all; clear; clc;
cd '/Users/ece.kaya/MATLAB/Githubbed/Misc/MTF';


kappa = .1;

maxtime = 60;
timeInc = .2;

nsamps = length(1:timeInc:maxtime);

stimmax = 1.5;
stimmin = 0;

% generate stim 
stim = [];
while 1
    
    stim = [stim repmat(stimmax,60,1)' repmat(stimmin,10,1)'];
    
    if length(stim) > nsamps
        stim = stim(1:nsamps);
        break
    end
end
figure;
stem(stim);

osc = table();
osc.pha_1 = stim';

figure;
set(gcf,'Position',[100 100 500 1000]);


subplot(3,1,1)

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
    osc.pha_1_t(k) = mod(t + pha_1, 2*pi);
    osc.pha_2_t(k) = mod(t + pha_2, 2*pi);
    
    osc.wave_1(k) = pha_1;
    osc.wave_2(k) = sin(t + pha_2);
    
    make_dots(pha_1); hold on;
    make_dots(pha_2); hold off;
    drawnow
    
    k = k+1;
end

subplot(3,1,2)
plot(osc.t, osc.pha_1, 'linewidth', 2); hold on;
plot(osc.t, osc.wave_2, 'linewidth', 2);
EK_plotlabels('time', 'wave (sin(t + pha_2))', '',14);

subplot(3,1,3)
plot(osc.t, osc.pha_1_t, 'linewidth', 2); hold on;
plot(osc.t, osc.pha_2_t, 'linewidth', 2);
EK_plotlabels('time', 'phases(t + pha_n)', '',14);
%% stim hilbert

close all; clear; clc;
cd '/Users/ece.kaya/MATLAB/Githubbed/Misc/MTF';


kappa = .1;

maxtime = 3;
timeInc = .01;

nsamps = length(1:timeInc:maxtime);

stimmax = 1;
stimmin = -1;

% generate stim 
stim = [];
while 1
    
    stim = [stim repmat(stimmax,30,1)' repmat(stimmin,10,1)'];
    
    if length(stim) > nsamps
        stim = stim(1:nsamps);
        break
    end
end

% figure;
% subplot(2,1,1)
% stem(stim);
% arcstim = asin(stim);
% subplot(2,1,2)
% plot(arcstim);

osc = table();
% osc.pha_1 = arcstim';

figure;
set(gcf,'Position',[100 100 500 1000]);


subplot(3,1,1)

k = 1;
for t = 1:timeInc:maxtime % sampling rate
    
    osc.t(k) = t;
    
    if k == 1
        pha_1 = 0;
        pha_2 = 0;
        
    else
%         pha_1 = osc.pha_1(k); % stim
        pha_1 = asin(stim(k));
        pha_2 = .2 + osc.pha_2(k-1) + kappa * sin(osc.pha_1(k-1) - osc.pha_2(k-1));
        
    end
    
    osc.pha_1(k) = pha_1;
    osc.pha_2(k) = pha_2;
%     osc.pha_1_t(k) = mod(t + pha_1, 2*pi);
%     osc.pha_2_t(k) = mod(t + pha_2, 2*pi);
    
    osc.wave_1(k) = sin(pha_1);
    osc.wave_2(k) = sin(t + pha_2);
    
    make_dots(pha_1); hold on;
    make_dots(pha_2); hold off;
    drawnow
    
    k = k+1;
end

subplot(3,1,2)
plot(osc.t, osc.wave_1, 'linewidth', 2); hold on;
plot(osc.t, osc.wave_2, 'linewidth', 2);
EK_plotlabels('time', 'wave (sin(t + pha_2))', '',14);

subplot(3,1,3)
plot(osc.t, osc.pha_1, 'linewidth', 2); hold on;
plot(osc.t, osc.pha_2, 'linewidth', 2);
EK_plotlabels('time', 'phases(t + pha_n)', '',14);

%% kappas ???
close all; clear; clc;

clrs = lines(10);


K = .01;

n = 5;      % number of osc
times = [];

waves = [];

kappas = [0    0    0    0    0 ;...      % 1st osc (stim) is independent
          K    0    K    K    K;...      % 2nd (auditory) 
          K    K    0    K    K;...
          K    K    K    0    K;...          
          K    K    K    K    0]; 

      
phases = [(1:n)'/n*2*pi]';

figure;
set(gcf,'Position',[100 100 500 1000])
subplot(2,1,1)
k = 1;
for t = 1:.2:20 % time
    
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
    plot(times, waves(:,i), 'Color', clrs(i,:), 'linewidth', 2); hold on;
end



%% kappas ???
close all; clear; clc;

clrs = lines(10);


K = .1;

n = 2;      % number of osc
times = [];

waves = [];

kappas = [0    K ;...      % 1st osc (stim) is independent
          K    0]; 

      
phases = [(1:n)'/n*2*pi]';

figure;
set(gcf,'Position',[100 100 500 1000])
subplot(2,1,1)
k = 1;
for t = 1:.2:20 % time
    
    for i = 1:n
        if k == 1
            pha = phases(k,i);
            
        else
            
            sum_sin_pha_diff = 0;
            for j = 1:n
                sum_sin_pha_diff = sum_sin_pha_diff + ...
                kappas(i,j) * sin(phases(k-1,j) - phases(k-1,i));
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
    plot(times, waves(:,i), 'Color', clrs(i,:), 'linewidth', 2); hold on;
end

%% basic & omegas (intrin freq) (OK) 

close all; clear; clc;

kappa = .01;

n = 3;      % number of osc
times = [];
phases = [];
waves = [];

for i = 1:n
    phases(1,i) = 2*pi*rand(1);
    omegas(i) = randn(1)*.1;
end


figure;
set(gcf,'Position',[100 100 500 1000])
subplot(2,2,1);
for i = 1:n
   plot(i, omegas(i), "Marker",".","MarkerSize",50); hold on;
end
EK_plotlabels('osc', 'omega','',14);

k = 1;
for t = 1:.2:50 % time
    
    for i = 1:n
        
        
        
        
        
        if k == 1
            pha = phases(k,i);
            
        else
            
            sum_sin_pha_diff = 0;
            for j = [1:i-1 i+1:n]
                sum_sin_pha_diff = sum_sin_pha_diff + ...
                    sin(phases(k-1,j) - phases(k-1,i));
            end
            pha = omegas(i) + phases(k-1,i) + kappa * sum_sin_pha_diff;
            
        end
        
        times(k) = t;
        phases(k,i) = pha;
        waves(k,i) = sin(t + pha);
        
        subplot(2,2,2)
        make_dots(t + pha); 
        drawnow
        hold on;
%         EK_plotlabels('coord', 'coord','',14);
    end
    
    
    hold off;
    k = k+1;
end

subplot(2,2,[3 4])
for i = 1:n
    plot(times, waves(:,i)); hold on;
end
EK_plotlabels('time', 'phases','',14);
%% basic & omegas (intrin freq) (OK) ++ kappa shutdown 

close all; clear; clc;

kappa = .03;
t_kappa_end = 20;
kappalist = [];

n = 3;      % number of osc
times = [];
phases = [];
waves = [];

for i = 1:n
    phases(1,i) = 2*pi*rand(1);
    omegas(i) = abs(randn(1)*.1);
%     omegas(i) = 0;
end


figure;
set(gcf,'Position',[100 100 800 1000])
subplot(2,2,1);
for i = 1:n
   stem(i, omegas(i), "Marker",".","MarkerSize",50, "linewidth", 3); hold on;
end
EK_plotlabels('osc', 'omega','',14);



k = 1;
for t = 1:.2:30 % time
    
    for i = 1:n
        
        if t > t_kappa_end 
            kappa = 0;
        end
        
        
        
        if k == 1
            pha = phases(k,i);
            
        else
            
            sum_sin_pha_diff = 0;
            for j = [1:i-1 i+1:n]
                sum_sin_pha_diff = sum_sin_pha_diff + ...
                    sin(phases(k-1,j) - phases(k-1,i));
            end
            pha = omegas(i) + phases(k-1,i) + kappa * sum_sin_pha_diff;
            
        end
        
        times(k) = t;
        phases(k,i) = pha;
        waves(k,i) = sin(t + pha);
        kappalist(k) = kappa;
        
        subplot(2,2,2)
        make_dots(t + pha); 
        drawnow
        hold on;
%         EK_plotlabels('coord', 'coord','',14);
    end
    
    
    hold off;
    k = k+1;
end

subplot(2,2,[3 4])
for i = 1:n
    plot(times, waves(:,i),'LineWidth',2); hold on;
end
plot(times,double(kappalist' > 0),'LineWidth',2);
EK_plotlabels('time', 'phases','',14);
%% basic

close all; clear; clc;

kappa = .05;

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
for t = 1:.2:20 % time
    
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


