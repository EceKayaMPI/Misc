


f = 5;
fs = 44100;  
Ts = 1/44100;              % sampling period
ncycs = f/dur;          % number of cycles

Tvec = 1/fs:Ts:dur;     % time vector. <<< changed to ones(..)



% mod_sig = cos(2*pi*fm*t+pi); 


fi  = f;
for ii = 1:ncycs
    cycdur = 1/fi(ii);
    cycT   = Ts:Ts:cycdur;
    test = [test cos(2*pi*fi(ii)*cycT+pi)];
end

% while length(test) ~= length(Tvec) %% keeps generating stim until done
%     fi  = f + (SD*f) .* randn(1,ncycs);
%     fid = sum(1./fi);
%     test = [];
%     for ii = 1:ncycs
%         cycdur = 1/fi(ii);
%         cycT   = Ts:Ts:cycdur;
%         test = [test cos(2*pi*fi(ii)*cycT+pi)];
%     end
% end
% figure, plot(test)

