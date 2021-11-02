
%%%%%%%%%%%%%%%%%%%%%%%%%
% this script is used to create Figure 1 in the comment by Obleser, Henry,
% & Lakatos, What do we talk about when we talk about rhythm;
% a formal comment on Breska A, Deouell LY (2017) Neural mechanisms of rhythm-based
% temporal prediction: Delta phase-locking reflects temporal predictability but not rhythmic entrainment.
% PLoS Biol 15(2): e2001665.
% doi:10.1371/journal.pbio.2001665
%
% 
% (c) PETER LAKATOS, plakatos0@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%

%%%% paradigms 1-4 and 6-7 are based on the paradigm of Breska and Deouell, 2017
paradigm = 	1; % 1 - regular short SOA, 2 - regular long SOA, 3 - 500-900 ms, short random SOA, 4 - 1100-1500 ms, long random SOA, 5 - truly random SOA example, 6 - short SOA paired stimulation, 7 - long SOA paired stimulation







%%% 0 to disable, 1 to enable specific ERP components
bool_p1n1   = 0;
bool_cnv    = 0;
bool_p3     = 0;

%%%% amplitude scaling for ongoing oscillation and ERP components, of course this influences the results
amp_ongoing         = 1;
amp_p1n1            = 0.5;
amp_p3              = 1;
amp_cnv             = 1.5;

%%%% ongoing oscillation parameters
segment_length      = 60;         %%% !!! the ratio of segment length and modulation speed has to be an integer
modulation_speed    = 2;
frq1                = 1.7;

sampling_rate       = 1000;

%%%% stimulus parameters

stim_count           = 7;
stim_start           = 20;
number_of_blocks    = 100;
soa                 = 1.5;
soa_rand1           = 0.5:0.05:0.9; %% breska fast "random"
soa_rand2           = 1.1:0.05:1.5; %% breska slow "random"
soa_rand3           = 0.5:0.05:2; %% very random
soa_rand4           = 1.5:0.05:1.9; %% breska fast pairs
soa_rand5           = 1.9:0.1:2.7; %% breska slow pairs
if paradigm ==0
    soa = soa;
elseif paradigm ==1
    soa                 = 0.7;
elseif paradigm ==2
    soa                 = 1.3;
elseif paradigm ==3
    soa  = mean(soa_rand1);
elseif  paradigm ==4
    soa  = mean(soa_rand2);
elseif  paradigm ==5
    soa  = mean(soa_rand3);
elseif  paradigm ==6
    soa  = 0.7;
    stim_count = stim_count*2;
elseif  paradigm ==7
    soa  = 1.3;
    stim_count = stim_count*2;
end


%%%% phase reset, entrainment parameters
targetphase         = pi/2;     %%% for phase reset
entrainment_memory  = 4*soa;    %%% time (in seconds) it takes to return to default frequency (frq1)
entrainment_ratio1  = 1;        %%% the amount with which the wavelength becomes similar to the SOA on every stimulus repetition. 1 is full entrainment
entrainment_ratio2  = 0.5;      %%% the amount with which the wavelength becomes similar to the SOA on every stimulus repetition. 1 is full entrainment
response_onset      = 0.01;     %%% poststimulus time at which reset occurs and ERP starts
reset_duration     = 0;     %%% duration of phase reset, set it to 0 for instantaneous reset

if reset_duration == 0
    reset_duration = 1/sampling_rate;
end

%%% epoching parameters
epoch_tframe_start      = [-10 20];
epoch_tframe_end        = [-20 10];
image_tframe            = [-10 10];

%%%% ERP components
% p1-n1-p2
t_p1 = 0.014:1/sampling_rate:0.04;
a_p1 = 0.5;
a_n1 = 1;
w_p1 = a_p1*cos(pi+2*pi*18.*t_p1);
t_n1 = 0.03:1/sampling_rate:0.156;
w_n1 = a_n1*cos(2*pi*8.*t_n1);
w_p1n1=amp_p1n1*[w_p1 w_n1];    %% p1n1
t_p1n1=(0:length(w_p1n1)-1)/sampling_rate;

% p3
t_p3 = 0.25:1/sampling_rate:0.416;
a_p3 = 1;
w_p3 = a_p3*cos(2*pi*3.*t_p3);
w_p3 = amp_p3*[zeros(1,length(0:1/sampling_rate:0.2)),w_p3];  %% p3
t_p31=(0:length(w_p3)-1)/sampling_rate;

% cnv
t_cnv1 = 0:1/sampling_rate:(soa-0.1);
t_cnv2 = 0:1/sampling_rate:0.147;
a_cnv = 1;
w_cnv1 = a_cnv*cos(pi/2+2*pi*(0.3/soa).*t_cnv1);
w_cnv2 = a_cnv*cos(pi/2+pi/2+2*pi*1.7.*t_cnv2);
w_cnv = amp_cnv*[zeros(1,length(w_p1n1)) w_cnv1 w_cnv2];  %% cnv
t_cnv=(0:length(w_cnv)-1)/sampling_rate;


stim_end        = stim_start+(stim_count-1)*soa;

t               = 0:1/sampling_rate:segment_length;

rand_phase_delay = [rand(1,ceil(number_of_blocks/2))*pi, rand(1,ceil(number_of_blocks/2))*-pi ];  %%% this is needed to avoid "synchronization" at oscillation initiation

alltrials = zeros(number_of_blocks,stim_count);
soa_rands=[];
for blockcik = 1:number_of_blocks
    if paradigm == 0 ||  paradigm == 1 || paradigm == 2
        alltrials(blockcik,:)=stim_start:soa:stim_end;
    end
    if paradigm == 3
        soa_rand = randsample(soa_rand1,stim_count,'true');
        soa_rands=[soa_rands soa_rand];
        %%% to randomize onset, use this option
        % alltrials(blockcik,1)=stim_start-soa+soa_rand(1);
        %%% otherwise
        alltrials(blockcik,1)=stim_start;
        for i1=2:stim_count
            alltrials(blockcik,i1)=alltrials(blockcik,i1-1)+soa_rand(i1);
        end
    end
    if paradigm == 4
        soa_rand = randsample(soa_rand2,stim_count,'true');
        soa_rands=[soa_rands soa_rand];
        %%% to randomize onset, use this option
        % alltrials(blockcik,1)=stim_start-soa+soa_rand(1);
        %%% otherwise
        alltrials(blockcik,1)=stim_start;
        for i1=2:stim_count
            alltrials(blockcik,i1)=alltrials(blockcik,i1-1)+soa_rand(i1);
        end
    end
    if paradigm == 5
        soa_rand = randsample(soa_rand3,stim_count,'true');
        soa_rands=[soa_rands soa_rand];
        %%% to randomize onset, use this option
        % alltrials(blockcik,1)=stim_start-soa+soa_rand(1);
        %%% otherwise
        alltrials(blockcik,1)=stim_start;
        for i1=2:stim_count
            alltrials(blockcik,i1)=alltrials(blockcik,i1-1)+soa_rand(i1);
        end
    end
    if paradigm == 6
        
        sr = randsample(soa_rand4,stim_count/2,'true');
        a=0;
        for i1=1:length(sr)
            a=a+1;
            soa_rand(a)=sr(i1);
            a=a+1;
            soa_rand(a)=0.7;
        end
        soa_rands=[soa_rands soa_rand];
        %%% to randomize onset, use this option
        % alltrials(blockcik,1)=stim_start-soa+soa_rand(1);
        %%% otherwise
        alltrials(blockcik,1)=stim_start;
        for i1=2:length(soa_rand)
            alltrials(blockcik,i1)=alltrials(blockcik,i1-1)+soa_rand(i1);
        end
    end
    if paradigm == 7
        
        sr = randsample(soa_rand5,stim_count/2,'true');
        a=0;
        for i1=1:length(sr)
            a=a+1;
            soa_rand(a)=sr(i1);
            a=a+1;
            soa_rand(a)=1.3;
        end
        soa_rands=[soa_rands soa_rand];
        %%% to randomize onset, use this option
        % alltrials(blockcik,1)=stim_start-soa+soa_rand(1);
        %%% otherwise
        alltrials(blockcik,1)=stim_start;
        for i1=2:length(soa_rand)
            alltrials(blockcik,i1)=alltrials(blockcik,i1-1)+soa_rand(i1);
        end
    end
end

trig1_pos=[];
for blockcik = 1:number_of_blocks
    for i1=1:stim_count
        trig1_pos(blockcik,i1)   = find(t<=alltrials(blockcik,i1), 1, 'last' );
    end
end

t_frqs              = zeros(number_of_blocks, length(t));
y_approxs           = zeros(number_of_blocks, length(t));
ph_origs            = zeros(number_of_blocks, length(t));
y_approx_resets     = zeros(number_of_blocks, length(t));
ph_resets           = zeros(number_of_blocks, length(t));
t_frq_entraineds    = zeros(number_of_blocks, length(t));
ph_entraineds       = zeros(number_of_blocks, length(t));
y_approx_entraineds = zeros(number_of_blocks, length(t));
t_frq_entrained50s  = zeros(number_of_blocks, length(t));
ph_entrained50s     = zeros(number_of_blocks, length(t));
y_approx_entrained50s = zeros(number_of_blocks, length(t));

for blockcik = 1:number_of_blocks
    
    trig1=alltrials(blockcik,:);
    
    %%% ongoing signal
    x               = frq1 + normrnd(0,0.2*frq1,segment_length/modulation_speed+1,1)'; %%% frequency changing randomly every modulation_speed second
    x_time          = 0:modulation_speed:segment_length;
    t_frq           = interp1(x_time,x,t,'pchip'); %%% instantaneous frequency
    y_wrong         = amp_ongoing*cos(2*pi*t_frq.*t);
    y_approx        = amp_ongoing*cos(2*pi*cumsum(t_frq)/sampling_rate+rand_phase_delay(blockcik));  %%% this is the base oscillation calculated the right way
    
    t_frqs(blockcik,:)=t_frq;
    y_approxs(blockcik,:)=y_approx;
    
    %%% phase reset signal
    ph                      = 2*pi*cumsum(t_frq)/sampling_rate+rand_phase_delay(blockcik);
    ph_origs(blockcik,:)    = ph - 2*pi*floor( (ph+pi)/(2*pi) );
    for i1=1:length(trig1)
        x1 = targetphase-ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+response_onset*sampling_rate+reset_duration*sampling_rate);
        ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)=ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)+x1;
    end
    ph      = ph - 2*pi*floor( (ph+pi)/(2*pi) );
    ph_resets(blockcik,:) = ph;
    y_approx_reset      = amp_ongoing*cos(ph_resets(blockcik,:));  %%% this is the base oscillation calculated the right way
    for i1=1:length(trig1)
        x1 = linspace(y_approx_reset(trig1_pos(blockcik,i1)+response_onset*sampling_rate),y_approx_reset(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate),reset_duration*sampling_rate+1);
        y_approx_reset(trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate)=x1;
    end
    
    t2 = t;
    y_approx_resets(blockcik,:) = y_approx_reset;
    
    
    %%% total entrainment signal
    z1                  = find(x_time<trig1(1));                                     % random frequency fluctuation before entrainment
    z2                  = find(x_time>trig1(end)+entrainment_memory);  % random frequency fluccutation after entrainment and entrainment memory
    
    
    xx=[];
    xx_time=[];
    a = 0;
    for i1=1:length(trig1)
        if i1==1
            a=a+1;
            last_ongoing_frq = x(z1(end));
            xx(a)=last_ongoing_frq;
            xx_time(a)=trig1(1);
        else
            a=a+1;
            xx(a)=xx(a-1);
            xx_time(a)=trig1(i1)+response_onset-1/sampling_rate;
            a=a+1;
            reprate             = 1/(trig1(i1)-trig1(i1-1));
            xx(a)=xx(a-1)+(reprate-xx(a-1))*entrainment_ratio1;
            xx_time(a)=xx_time(a-1)+1/sampling_rate;
        end
    end
    a=a+1;
    xx(a)=xx(a-1)+(reprate-xx(a-1))*entrainment_ratio1;
    xx_time(a)=trig1(end)+soa*1+response_onset;
    
    x_entrained         = [x(z1) xx x(z2)];
    x_entrained_time     = [x_time(z1) xx_time x_time(z2)];
    t_frq_entrained     = interp1(x_entrained_time,x_entrained,t,'pchip'); %%% instantaneous frequency
    t_frq_entraineds(blockcik,:)=t_frq_entrained;
    
    %%%% have to also reset the frequency matched oscillation
    ph      = 2*pi*cumsum(t_frq_entrained)/sampling_rate;
    for i1=1:length(trig1)
        x1 = targetphase-ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+response_onset*sampling_rate+reset_duration*sampling_rate);
        ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)=ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)+x1;
    end
    ph      = ph - 2*pi*floor( (ph+pi)/(2*pi) );
    ph_entraineds(blockcik,:) = ph;
    y_approx_entrained      = amp_ongoing*cos(ph_entraineds(blockcik,:));  %%% this is the base oscillation calculated the right way
    for i1=1:length(trig1)
        x1 = linspace(y_approx_entrained(trig1_pos(blockcik,i1)+response_onset*sampling_rate),y_approx_entrained(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate),reset_duration*sampling_rate+1);
        y_approx_entrained(trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate)=x1;
    end
    
    % correction: until 2nd stimulus just use phase reset signal
    z=find(t<trig1(2)+response_onset);
    y_approx_entrained(z) = y_approx_reset(z);
    
    y_approx_entraineds(blockcik,:)=y_approx_entrained;
    
    %%% 50% entrainment signal, only entrainment ratio is changed
    
    z1                  = find(x_time<trig1(1));                                      % random frequency fluctuation before entrainment
    z2                  = find(x_time>trig1(end)+entrainment_memory);  % random frequency fluccutation after entrainment and entrainment memory
    
    
    xx=[];
    xx_time=[];
    a = 0;
    for i1=1:length(trig1)
        if i1==1
            a=a+1;
            last_ongoing_frq = x(z1(end));
            xx(a)=last_ongoing_frq;
            xx_time(a)=trig1(1);
        else
            
            a=a+1;
            xx(a)=xx(a-1);
            xx_time(a)=trig1(i1)-1/sampling_rate+response_onset;
            a=a+1;
            reprate             = 1/(trig1(i1)-trig1(i1-1));
            xx(a)=xx(a-1)+(reprate-xx(a-1))*entrainment_ratio2;
            xx_time(a)=trig1(i1)+response_onset;
        end
    end
    a=a+1;
    xx(a)=xx(a-1)+(reprate-xx(a-1))*entrainment_ratio2;
    xx_time(a)=trig1(end)+soa*1+response_onset;
    
    x_entrained50         = [x(z1) xx x(z2)];
    x_entrained50_time     = [x_time(z1) xx_time x_time(z2)];
    t_frq_entrained50     = interp1(x_entrained50_time,x_entrained50,t,'pchip'); %%% instantaneous frequency
    t_frq_entrained50s(blockcik,:)=t_frq_entrained50;
    
    %%%% have to also reset the frequency matched oscillation
    ph      = 2*pi*cumsum(t_frq_entrained50)/sampling_rate;
    for i1=1:length(trig1)
        x1 = targetphase-ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+response_onset*sampling_rate+reset_duration*sampling_rate);
        ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)=ph(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate:end)+x1;
    end
    ph      = ph - 2*pi*floor( (ph+pi)/(2*pi) );
    ph_entrained50s(blockcik,:) = ph;
    y_approx_entrained50      = amp_ongoing*cos(ph_entrained50s(blockcik,:));  %%% this is the base oscillation calculated the right way
    
    for i1=1:length(trig1)
        x1 = linspace(y_approx_entrained50(trig1_pos(blockcik,i1)+response_onset*sampling_rate),y_approx_entrained50(trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate),reset_duration*sampling_rate+1);
        y_approx_entrained50(trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+reset_duration*sampling_rate)=x1;
    end
    
    % correction: until 2nd stimulus just use phase reset signal
    z=find(t<trig1(2)+response_onset);
    y_approx_entrained50(z) = y_approx_reset(z);
    
    y_approx_entrained50s(blockcik,:)=y_approx_entrained50;
    
    
    
end
if bool_p1n1 == 1
    for blockcik = 1:number_of_blocks
        for i1=1:length(trig1)
            y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)=y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)+w_p1n1;
            y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)=y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)+w_p1n1;
            y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)=y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)+w_p1n1;
            y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)=y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p1n1)-1)+w_p1n1;
        end
    end
end

if bool_p3 == 1
    for blockcik = 1:number_of_blocks
        i1=length(trig1);
        y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)=y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)+w_p3;
        y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)=y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)+w_p3;
        y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)=y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)+w_p3;
        y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)=y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_p3)-1)+w_p3;
    end
end

if bool_cnv == 1
    for blockcik = 1:number_of_blocks
        i1=length(trig1)-1;
        y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)=y_approxs(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)+w_cnv;
        y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)=y_approx_resets(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)+w_cnv;
        y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)=y_approx_entraineds(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)+w_cnv;
        y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)=y_approx_entrained50s(blockcik,trig1_pos(blockcik,i1)+response_onset*sampling_rate:trig1_pos(blockcik,i1)+response_onset*sampling_rate+length(w_cnv)-1)+w_cnv;
    end
end

%%%% Breska phase estimation, Breska used n = 4 (24 dB/oct), but I got an error

y_approxs_f = zeros(size(y_approxs));
y_approx_resets_f     = zeros(size(y_approxs));
y_approx_entraineds_f   = zeros(size(y_approxs));
y_approx_entrained50s_f = zeros(size(y_approxs));

y_approxs_f_ph           = zeros(size(y_approxs));
y_approx_resets_f_ph     = zeros(size(y_approxs));
y_approx_entraineds_f_ph   = zeros(size(y_approxs));
y_approx_entrained50s_f_ph = zeros(size(y_approxs));

n = 3;
Wn = [0.5 3]/(sampling_rate/2);
[b,a] = butter(n,Wn);
for blockcik = 1:number_of_blocks
    for i1=1:length(trig1)
        y_approxs_f(blockcik,:)             = filtfilt(b,a,y_approxs(blockcik,:));
        y_approx_resets_f(blockcik,:)       = filtfilt(b,a,y_approx_resets(blockcik,:));
        y_approx_entraineds_f(blockcik,:)   = filtfilt(b,a,y_approx_entraineds(blockcik,:));
        y_approx_entrained50s_f(blockcik,:) = filtfilt(b,a,y_approx_entrained50s(blockcik,:));
        
        y_approxs_f_ph(blockcik,:)             = angle(hilbert(y_approxs_f(blockcik,:)));
        y_approx_resets_f_ph(blockcik,:)       = angle(hilbert(y_approx_resets_f(blockcik,:)));
        y_approx_entraineds_f_ph(blockcik,:)   = angle(hilbert(y_approx_entraineds_f(blockcik,:)));
        y_approx_entrained50s_f_ph(blockcik,:) = angle(hilbert(y_approx_entrained50s_f(blockcik,:)));
        
    end
    
end


%%%%%%%%%% epoching

dt = 1/sampling_rate;
trig01 = trig1_pos(:,1);
trig02 = trig1_pos(:,end);
x1s = round(epoch_tframe_start(1)/dt);
x2s = round(epoch_tframe_start(2)/dt);
x1e = round(epoch_tframe_end(1)/dt);
x2e = round(epoch_tframe_end(2)/dt);
eeg_time_start = epoch_tframe_start(1):dt:epoch_tframe_start(2);
eeg_time_end = epoch_tframe_end(1):dt:epoch_tframe_end(2);
eeg_start_ongoing   = zeros(blockcik,length(eeg_time_start));
eeg_end_ongoing     = zeros(blockcik,length(eeg_time_end));
eeg_start_reset   = zeros(blockcik,length(eeg_time_start));
eeg_end_reset     = zeros(blockcik,length(eeg_time_end));
eeg_start_entrained   = zeros(blockcik,length(eeg_time_start));
eeg_end_entrained     = zeros(blockcik,length(eeg_time_end));
eeg_start_entrained50   = zeros(blockcik,length(eeg_time_start));
eeg_end_entrained50     = zeros(blockcik,length(eeg_time_end));
for i1=1:length(trig01)
    eeg_start_ongoing(i1,:)     = y_approxs(i1,trig01(i1)+x1s:trig01(i1)+x2s);
    eeg_end_ongoing(i1,:)       = y_approxs(i1,trig02(i1)+x1e:trig02(i1)+x2e);
    eeg_start_reset(i1,:)     = y_approx_resets(i1,trig01(i1)+x1s:trig01(i1)+x2s);
    eeg_end_reset(i1,:)       = y_approx_resets(i1,trig02(i1)+x1e:trig02(i1)+x2e);
    eeg_start_entrained(i1,:)     = y_approx_entraineds(i1,trig01(i1)+x1s:trig01(i1)+x2s);
    eeg_end_entrained(i1,:)       = y_approx_entraineds(i1,trig02(i1)+x1e:trig02(i1)+x2e);
    eeg_start_entrained50(i1,:)     = y_approx_entrained50s(i1,trig01(i1)+x1s:trig01(i1)+x2s);
    eeg_end_entrained50(i1,:)       = y_approx_entrained50s(i1,trig02(i1)+x1e:trig02(i1)+x2e);
end


%%%%%%%%% ITC calculation based on real and measured phases

%%%%% real
for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=ph_origs(blockcik,trig1_pos(blockcik,i1));
    end
    [p_origs(1,i1), r_origs(1,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=ph_resets(blockcik,trig1_pos(blockcik,i1));
    end
    [p_resets(1,i1), r_resets(1,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=ph_entraineds(blockcik,trig1_pos(blockcik,i1));
    end
    [p_entraineds(1,i1), r_entraineds(1,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=ph_entrained50s(blockcik,trig1_pos(blockcik,i1));
    end
    [p_entrained50s(1,i1), r_entrained50s(1,i1)] = rayleigh(phs');
end

%%%%% measured
for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=y_approxs_f_ph(blockcik,trig1_pos(blockcik,i1));
    end
    [p_origs(2,i1), r_origs(2,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=y_approx_resets_f_ph(blockcik,trig1_pos(blockcik,i1));
    end
    [p_resets(2,i1), r_resets(2,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=y_approx_entraineds_f_ph(blockcik,trig1_pos(blockcik,i1));
    end
    [p_entraineds(2,i1), r_entraineds(2,i1)] = rayleigh(phs');
end

for i1=1:length(trig1)
    phs=[];
    for blockcik = 1:number_of_blocks
        phs(blockcik)=y_approx_entrained50s_f_ph(blockcik,trig1_pos(blockcik,i1));
    end
    [p_entrained50s(2,i1), r_entrained50s(2,i1)] = rayleigh(phs');
end



%%% plotting
fsize = 8;
seltrials = 1;



curfig = figure;
set(curfig,'position',[100  50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto','name',['entrainment simulator - fig1'])
axes('Position',[0 0.985 1 0.2],'Visible','off');
text(0.5,0,['entrainment simulator - SOA ' num2str(soa*1000) ' ms, fig1, ERPs'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')

subplot(5,3,1)
plot(t,t_frqs(seltrials,:),t,t_frq_entraineds(seltrials,:),t,t_frq_entrained50s(seltrials,:));
set(gca,'xlim',[0 segment_length],'fontsize',fsize)
hold on
for i1=1:length(trig1)
    line([mean(alltrials(seltrials,i1),1) mean(alltrials(seltrials,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
yl(1,:)=get(gca,'ylim');
legend('ongoing','100% entrained', '50% entrained')
title('frequency')

subplot(5,3,4)
plot(t,y_approxs(seltrials,:));
hold on
for i1=1:length(trig1)
    line([mean(alltrials(seltrials,i1),1) mean(alltrials(seltrials,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
set(gca,'xlim',[0 segment_length],'fontsize',fsize)
yl(2,:)=get(gca,'ylim');
title ('ongoing, single trial')

subplot(5,3,7)
plot(t2,y_approx_resets(seltrials,:));
hold on
for i1=1:length(trig1)
    line([mean(alltrials(seltrials,i1),1) mean(alltrials(seltrials,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
set(gca,'xlim',[0 segment_length],'fontsize',fsize)
yl(3,:)=get(gca,'ylim');
title ('phase reset only, single trial')

subplot(5,3,10)
plot(t,y_approx_entraineds(seltrials,:));
hold on
for i1=1:length(trig1)
    line([mean(alltrials(seltrials,i1),1) mean(alltrials(seltrials,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
set(gca,'xlim',[0 segment_length],'fontsize',fsize)
yl(4,:)=get(gca,'ylim');
title ('100% wavelength modulation, single trial')

subplot(5,3,13)
plot(t,y_approx_entrained50s(seltrials,:));
hold on
for i1=1:length(trig1)
    line([mean(alltrials(seltrials,i1),1) mean(alltrials(seltrials,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
set(gca,'xlim',[0 segment_length],'fontsize',fsize)
yl(5,:)=get(gca,'ylim');
title ('50% wavelength modulation, single trial')

subplot(5,3,2)
plot(t,mean(t_frqs,1),t,mean(t_frq_entraineds,1),t,mean(t_frq_entrained50s,1));
set(gca,'xlim',[mean(alltrials(:,1),1)+eeg_time_start(1) mean(alltrials(:,1),1)+eeg_time_start(end)],'fontsize',fsize)
set(gca,'ylim',yl(1,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1) mean(alltrials(:,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
legend('ongoing','100% entrained', '50% entrained')
title('mean frequency')

subplot(5,3,5)
%plot(t,mean(y_approxs,1));
plot(eeg_time_start,mean(eeg_start_ongoing(:,:),1));
set(gca,'xlim',[eeg_time_start(1) eeg_time_start(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,1),1) mean(alltrials(:,i1),1)-mean(alltrials(:,1),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('ongoing, average, start synched')

subplot(5,3,8)
%plot(t2,mean(y_approx_resets,1));
plot(eeg_time_start,mean(eeg_start_reset(:,:),1));
set(gca,'xlim',[eeg_time_start(1) eeg_time_start(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,1),1) mean(alltrials(:,i1),1)-mean(alltrials(:,1),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('phase reset only, average, start synched')


subplot(5,3,11)
%plot(t,mean(y_approx_entraineds,1));
plot(eeg_time_start,mean(eeg_start_entrained(:,:),1));
set(gca,'xlim',[eeg_time_start(1) eeg_time_start(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,1),1) mean(alltrials(:,i1),1)-mean(alltrials(:,1),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('100% wavelength modulation, average, start synched')

subplot(5,3,14)
% plot(t,mean(y_approx_entrained50s,1));
plot(eeg_time_start,mean(eeg_start_entrained50(:,:),1));
set(gca,'xlim',[eeg_time_start(1) eeg_time_start(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,1),1) mean(alltrials(:,i1),1)-mean(alltrials(:,1),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('50% wavelength modulation, average, start synched')

subplot(5,3,3)
plot(t,mean(t_frqs,1),t,mean(t_frq_entraineds,1),t,mean(t_frq_entrained50s,1));
set(gca,'xlim',[mean(alltrials(:,end),1)+eeg_time_end(1) mean(alltrials(:,end),1)+eeg_time_end(end)],'fontsize',fsize)
set(gca,'ylim',yl(1,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1) mean(alltrials(:,i1),1)],get(gca,'ylim'),'color',[1 0 0])
end
legend('ongoing','100% entrained', '50% entrained')
title('mean frequency')

subplot(5,3,6)
%plot(t,mean(y_approxs,1));
plot(eeg_time_end,mean(eeg_end_ongoing(:,:),1));
set(gca,'xlim',[eeg_time_end(1) eeg_time_end(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,end),1) mean(alltrials(:,i1),1)-mean(alltrials(:,end),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('ongoing, average, end synched')

subplot(5,3,9)
%plot(t,mean(y_approxs,1));
plot(eeg_time_end,mean(eeg_end_reset(:,:),1));
set(gca,'xlim',[eeg_time_end(1) eeg_time_end(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,end),1) mean(alltrials(:,i1),1)-mean(alltrials(:,end),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('reset only, average, end synched')

subplot(5,3,12)
%plot(t,mean(y_approxs,1));
plot(eeg_time_end,mean(eeg_end_entrained(:,:),1));
set(gca,'xlim',[eeg_time_end(1) eeg_time_end(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,end),1) mean(alltrials(:,i1),1)-mean(alltrials(:,end),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('100% wavelength modulation, average, end synched')

subplot(5,3,15)
%plot(t,mean(y_approxs,1));
plot(eeg_time_end,mean(eeg_end_entrained50(:,:),1));
set(gca,'xlim',[eeg_time_end(1) eeg_time_end(end)],'fontsize',fsize)
set(gca,'ylim',yl(2,:))
hold on
for i1=1:length(trig1)
    line([mean(alltrials(:,i1),1)-mean(alltrials(:,end),1) mean(alltrials(:,i1),1)-mean(alltrials(:,end),1)],get(gca,'ylim'),'color',[1 0 0])
end
title ('50% wavelength modulation, average, end synched')

curfig = figure;
set(curfig,'position',[100  50   1500   900],'color',[1 1 1],'InvertHardcopy','off','PaperPositionMode','auto','name',['entrainment simulator - fig3'])
axes('Position',[0 0.985 1 0.2],'Visible','off');
text(0.5,0,['entrainment simulator - SOA ' num2str(soa*1000) ' ms, fig2, ITC values'],'FontSize',11,'color',[0 0 0],'HorizontalAlignment','center','interpreter','none')


subplot(1,2,1)
plot(1:length(trig1),r_origs(1,:),1:length(trig1),r_resets(1,:),1:length(trig1),r_entraineds(1,:),1:length(trig1),r_entrained50s(1,:))
legend('ongoing','reset only','100% entrained', '50% entrained')
set(gca,'xlim',[0.5 length(trig1)+0.5])
title('ITC based on real phases')
ylabel('ITC')
xlabel('stimulus')

subplot(1,2,2)
plot(1:length(trig1),r_origs(2,:),1:length(trig1),r_resets(2,:),1:length(trig1),r_entraineds(2,:),1:length(trig1),r_entrained50s(2,:))
legend('ongoing','reset only','100% entrained', '50% entrained')
set(gca,'xlim',[0.5 length(trig1)+0.5])
title('ITC based on measured phases')
ylabel('ITC')
xlabel('stimulus')
