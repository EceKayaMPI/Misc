% 
% x = 1:5;
% 
% EK_dotplot(x)
% 

function EK_dotplot(x)
scatter(x, ones(1,length(x)),'LineWidth',1);
% set(gca,'LineWidth',0.75);
end