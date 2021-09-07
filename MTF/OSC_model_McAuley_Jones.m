
function [model_con, hit_list] = OSC_model_McAuley_Jones (condit_list,Wphi,Wp)

model_con = struct();

% Wphi = .5;            % phase correction weight 
                        % [0 = no phase correction, 1 = full phase reset]

% Wp =  .5;             % period correction weight
                        % [0 = no period corr. 1 = complete period corr
                        % where P matches curr IOI]

pha = 0;              % phase corr term (model's default phase)
per = 350;            % period corr term (model's default period)

C = 0;                % temporal contrast (to be estimated)


devdir_list = condit_list(:,end);
stim_IOIs_list = condit_list(:,1:end-1);

ncondits = size(stim_IOIs_list,1);

for c = 1:ncondits
    
    stim_IOIs = stim_IOIs_list(c,:);
    devdir = devdir_list(c);
    
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
    model_con(c).devdir = devdir;
    
    model_con(c).C_fin = C;
    
    if devdir > 0 
        if C > 0 
            model_con(c).hit = 1;
        else
            model_con(c).hit = 0;
        end
    else
        if C > 0 
            model_con(c).hit = 0;
        else
            model_con(c).hit = 1;
        end 
    end

end

hit_list = [model_con.hit].';

end