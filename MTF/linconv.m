
% clear, clc;
% 
% ar = [.10 .40 -.40 -.10];

function [longer, shorter] = linconv (ar)
for i = 1:length(ar)
x = ar(i);
    
    if x > 0 
        
        if x >= .25   
           longer(i) = x - .25;
           shorter(i) = .25 + (.50 - x);
           
        elseif x < .25
           longer(i) = .25 - x;
           shorter(i) = .25 + x;
        end
        
    else % x < 0
        
        if abs(x) >= .25
            longer(i) = .25 + (.50 + x);
            shorter(i) = abs(x) - .25;
            
        elseif abs(x) < .25
            longer(i) = abs(x) + .25;
            shorter(i) = .25 + x;
        end
    end
    
end


end 
