%% calculates by cycle
clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;  
f = 8;

%--------multi freq
testm = [];
% flist = randi(2*f,5,1);
flist = repmat(10,5,1);

ioilist = [100, 100, 100, 200, 100];

for ioi = 1:5

    cycdur = ioilist(ioi)/fs;
    freq = 1/cycdur;
    cycT   = Ts:Ts:cycdur;
    testm = [testm cos(2*pi*freq*cycT)];
end


% subplot(2,1,2)
plot(testm, 'linewidth', 2);

%% calculates by time (Ts)
clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;

wave = [];
ioilist = [100, 100, 100, 200, 100];
testm = [];

k = 1;
for ioi = 1:5

    cycdur = ioilist(ioi)/fs;
    freq = 1/cycdur;
    
    % calculates by adding full waves, calc by multiplying time vector 
    cycT   = Ts:Ts:cycdur;
    testm = [testm cos(2*pi*freq*cycT+1)];
    

    for t = Ts:Ts:cycdur % just 1 interval 
        
        % calculates by adding wave points, calc by multiplying time value
        % (t)
        wave(k) = cos(t*(2*pi*freq)+1);
        k = k+1;
    end
end

figure;
plot(wave); hold on;
plot(testm);

%% full frequency (period) correction 
clear all;clc;
close all;

fs = 1000;  
Ts = 1/fs;

wave = [];

ioilist = [100, 100, 100, 100, 100];
stim = zeros(200,1)';
for i = 1:length(ioilist)
    x = zeros(ioilist(i),1)';
    x(1) = 1;
    stim = [stim x];
end

figure; stem(stim); hold on;
% make stim VECTOR 


k = 1;
t = Ts;

freq = 5; % fs/samples << initial freq 

tbl = table();
a = 1;

while 1     
    
    if stim(k) == 1 % if you see a stimulus (sound)

        if a == 1
            ioi = k;
        else
            ioi = k - tbl.k(a-1);
        end

        
        freqnew = fs/ioi;
        disp(['freq was ' num2str(freq) ' --> ' num2str(freqnew)]);
        freq = freqnew;        
        
        contin_ok = input('continue? (1/0): ');
        if contin_ok ==0
            break;
        end
        
        tbl.ioi(a) = ioi;
        tbl.freq(a) = freq;
        tbl.k(a) = k;
        a = a+1;
    end
    
    
    
    % calculates by adding wave points, calc by multiplying time value
    % (t)
    wave(k) = cos(t*(2*pi*freq)+0);
    plot(wave);
    
    
    t = t+Ts;
    if k == length(stim)
        break
    end
    
    k = k+1;
end



%%
clear all;clc;
close all;


fs = 1000;  
Ts = 1/fs;   

dur = 1;
f = 2;
ncycs = dur*f;

%--------single freq
test = [];
fi  = f;
for ioi = 1:1%ncycs
    cycdur = 1/fi;
    cycT   = Ts:Ts:cycdur;
    test = [test cos(2*pi*fi*cycT)];
end


figure;

% subplot(2,1,1); 
hold on
plot(test, 'linewidth', 2);
% 
% %--------multi freq
% testm = [];
% flist = randi(2*f,ncycs,1);
% for ii = 1:ncycs
%     cycdur = 1/flist(ii);
%     cycT   = Ts:Ts:cycdur;
%     testm = [testm cos(2*pi*flist(ii)*cycT-1)];
% end

% 
% % subplot(2,1,2)
% plot(testm, 'linewidth', 2);


%%

figure;
omega = 0;
k = 1;
phases(1) = 0;
waves = [];

subplot(2,1,1)
for t = 0:.1:2*pi % time
    
    if k == 1
        pha = 0;
    else
    
    pha = omega + phases(k-1);
    
    end
    

    phases(k) = pha;
    waves(k) = cos(t + pha);
    
    make_dots(t + pha);
    drawnow
    
  
    k = k+1;
end

subplot(2,1,2)
plot (waves);

%%
function make_dots(theta, clr)
if nargin<2
    clr = [0 0 0];
end
plot(cos(theta), sin(theta), "Marker",".","MarkerSize",50, "MarkerFaceColor", clr);
set(gca,'xlim', [-1 1], 'ylim', [-1 1]);
end
