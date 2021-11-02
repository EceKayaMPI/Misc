
function [pha3] = EK_phase_matcher (f1, f2, pha1, pha2, t)
% [pha3] = EK_phase_matcher (5, 6, 3, 5, .5) 
% new wave should have same freq as f2!! 

% radpha1 = t*2*pi*f1+pha1;
radpha2 = t*2*pi*f2+pha2;

radphadif = t*2*pi*(f2-f1)+(pha2-pha1);
pha2match = radpha2-radphadif;

pha3 = pha2match/360;

end
%% example
% % clear;clc;close all
% % 
% % f1 = 5;
% % f2 = 6;
% % 
% % pha1 = 3;
% % pha2 = 5;
% % 
% % Ts = 1/1000;
% % tvec = Ts:Ts:1;
% % 
% % figure;
% % subplot(2,1,1)
% % plot(tvec, cos(tvec*2*pi*f1+pha1))
% % hold on;
% % plot(tvec, cos(tvec*2*pi*f2+pha2)) 
% % 
% % legend(...
% %     ['freq = ' num2str(f1) ', phase = ' num2str(pha1)],...
% %     ['freq = ' num2str(f2) ', phase = ' num2str(pha2)]);
% % 
% % subplot(2,1,2)
% % plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f1+pha1))
% % hold on;
% % plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f2+pha2)) 
% % 
% % radpha1 = tvec(end-500)*2*pi*f1+pha1;
% % radpha2 = tvec(end-500)*2*pi*f2+pha2;
% % 
% % radphadif = tvec(end-500)*2*pi*(f2-f1)+(pha2-pha1);
% % disp(radphadif);
% % disp(radpha2 - radpha1);
% % 
% % pha2match = radpha2-radphadif;
% % 
% % f3 = f2;
% % pha3 = pha2match/360;
% % 
% % plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f3+pha2match)) 
% % 
% % legend(...
% %     ['freq = ' num2str(f1) ', phase = ' num2str(pha1)],...
% %     ['freq = ' num2str(f2) ', phase = ' num2str(pha2)],...
% %     ['freq = ' num2str(f3) ', phase = ' num2str(pha3)]);