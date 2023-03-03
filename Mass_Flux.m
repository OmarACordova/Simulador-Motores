function m=Pressure_Reg(P_high,P_low, T_high,gamma, molar_mass)
    Critical_Press_Ratio=(.5+.5*gamma)^(-gamma/(gamma-1));
    R=8.314;
    if(P_low/P_high<=Critical_Press_Ratio)
        Press_ratio=Critical_Press_Ratio;
        m=P_high*((2*gamma/((gamma-1)*R*T_high))^.5)*( (Press_ratio)^(2/gamma)-(Press_ratio)^((gamma+1)/gamma))
    else
        Press_ratio=P_low/P_high;
        m=P_high*((2*gamma/((gamma-1)*R*T_high))^.5)*( (Press_ratio)^(2/gamma)-(Press_ratio)^((gamma+1)/gamma))
    end
end