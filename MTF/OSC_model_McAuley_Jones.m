
function [model_con] = OSC_model_McAuley_Jones (stim_IOIs_list,Wphi,Wp)

model_con = struct();

% Wphi = .5;            % phase correction weight 
                        % [0 = no phase correction, 1 = full phase reset]

% Wp =  .5;             % period correction weight
                        % [0 = no period corr. 1 = complete period corr
                        % where P matches curr IOI]

pha = 0;              % phase corr term (model's default phase)
per = 350;            % period corr term (model's default period)

C = 0;                % temporal contrast (to be estimated)

ncondits = size(stim_IOIs_list,1);

for c = 1:ncondits
    
    stim_IOIs = stim_IOIs_list(c,:);
    
     for i = 1:length(stim_IOIs)
        
        if mod((pha(i) + stim_IOIs(i)/per(i)),1) > .5
            C = mod((pha(i) + stim_IOIs(i)/per(i)),1) - 1;
        else
            C = mod((pha(i) + stim_IOIs(i)/per(i)),1);
        end
        
        pha(i+1) = (1 - Wphi)*C;
        per(i+1) = (1 + Wp*C)*per(i);
        
     end
    
    model_con(c).stim_IOIs = stim_IOIs;
    model_con(c).C_fin = C;
    

end

% C_fin_list = 

end