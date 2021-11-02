
clear;clc;close all

f1 = 7;
f2 = 3;

theta1 = 4;
theta2 = 5;

Ts = 1/1000;
tvec = Ts:Ts:1;

figure;
subplot(2,1,1)
plot(tvec, cos(tvec*2*pi*f1+theta1))
hold on;
plot(tvec, cos(tvec*2*pi*f2+theta2)) 

legend(...
    ['freq = ' num2str(f1) ', phase = ' num2str(theta1)],...
    ['freq = ' num2str(f2) ', phase = ' num2str(theta2)]);

subplot(2,1,2)
plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f1+theta1))
hold on;
plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f2+theta2)) 
title('waves after t = .5');

radpha1 = tvec(end-500)*2*pi*f1+theta1;
% radpha1 = 18.708;
radpha2 = tvec(end-500)*2*pi*f2+theta2;

radphadif = tvec(end-500)*2*pi*(f2-f1)+(theta2-theta1);
% radphadif = radpha2 - radpha1;
disp(radphadif);
disp(radpha2 - radpha1);

pha2match = radpha2-radphadif;

f3 = f2;
theta3 = pha2match/360;

plot(tvec(end-500:end), cos(tvec(end-500:end)*2*pi*f3+pha2match)) 

legend(...
    ['freq = ' num2str(f1) ', phase = ' num2str(theta1)],...
    ['freq = ' num2str(f2) ', phase = ' num2str(theta2)],...
    ['freq = ' num2str(f3) ', phase = ' num2str(theta3)]);