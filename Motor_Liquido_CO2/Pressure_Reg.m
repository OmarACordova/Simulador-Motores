function [dm_flux, T_low ] =Pressure_Reg(P_high,P_obj, P_low, T_high,gamma, R)
    Critical_Press_Ratio=(.5+.5*gamma)^(-gamma/(gamma-1));
    
    if P_high+30^5>P_obj && P_obj>P_low
        
        if(P_low/P_obj<=Critical_Press_Ratio)
            Press_ratio=Critical_Press_Ratio;
            dm_flux=P_high*((2*gamma/((gamma-1)*R*T_high))^.5)*( (Press_ratio)^(2/gamma)-(Press_ratio)^((gamma+1)/gamma));

        else
            Press_ratio=P_low/P_obj;
            dm_flux=P_high*((2*gamma/((gamma-1)*R*T_high))^.5)*( (Press_ratio)^(2/gamma)-(Press_ratio)^((gamma+1)/gamma));
        end

        T_low=T_high*subplus(real((P_low/P_high)^(1-1/gamma)));

    else
        dm_flux=0;
        T_low=T_high;
    end
end