
clear;clc;close all

f2 = 6;
theta2 = 0;

Ts = 1/1000;
tvec = Ts:Ts:1;

figure;

hold on;
plot(tvec, cos(tvec*2*pi*f2+theta2)) 

radpha1 = 28.708;
radpha2 = tvec(end-500)*2*pi*f2+theta2;

radphadif = radpha2 - radpha1;


pha2match = radpha2-radphadif;

f3 = f2;
theta3 = pha2match/360;

plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f3+theta3)) 


legend(['freq = ' num2str(f2) ', phase = ' num2str(theta2)],...
    ['freq = ' num2str(f3) ', phase = ' num2str(theta3)]);