

% figure
% for i = -360:360
%     a = i/360*pi;
%     plot(cos(a), sin(a), 'ok'); hold on;
%     
% end
% 
% 
% figure
% for i = -pi:.1:pi
%     a = i;
%     plot(cos(a), sin(a), 'ok'); hold on;
%     
% end
% 
% figure
% for i = -pi: .1 :pi
%     a = i-pi;
%     plot(cos(a), sin(a), 'ok'); hold on;
%     
% end

figure
for i = 0:.1:1000
    i = mod(i, 2*pi);
    a = i+pi;
    plot(cos(a), sin(a), 'ok'); hold on;
    
end

%%

