

function EK_quickplot (xaxis, yaxis, xString, yString, tits) 
% figure;
plot(xaxis,yaxis,'LineWidth',1);
set(gca,'FontSize',14);
set(gca,'LineWidth',0.75);
xlabel(xString);
ylabel(yString);
title(tits);
end