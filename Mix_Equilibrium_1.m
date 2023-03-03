function [T_mix, Ro_1, Ro_2]=Mix_Equilibrium_1(H, m1,m2,Vol,sust1,sust2,Ro1_seed, T_seed)
    T_guess=T_seed;
    while T_guess<T_seed*1.1
        T_guess=T_guess+.1;
        Ro1_guess=Ro1_seed;
        while Ro1_guess*1.00001<1.1*Ro1_seed
            Ro1_guess=Ro1_guess+Ro1_seed*.0001;
            %[Ro1_seed Ro1_guess m1 Vol (m1/Ro1_guess)]    
            Ro2_guess=m2/(Vol-(m1/Ro1_guess));
            h1=py.CoolProp.CoolProp.PropsSI('H','T',T_guess,'D',Ro1_guess, sust1);
            h2=py.CoolProp.CoolProp.PropsSI('H','T',T_guess,'D',Ro2_guess, sust2);
            H_guess=h1*m1+h2*m2;
            %[Ro1_seed Ro1_guess Ro2_guess T_seed T_guess ]    

            if(H_guess/H < 1.001 && H_guess/H > 0.999)
                T_mix=T_guess;
                Ro_1=Ro1_guess;
                Ro_2=Ro2_guess;
                [Ro1_seed Ro_1 Ro_2 T_seed T_guess ]    
                a='Listo'
                %break
            end
        end
        
        if(H_guess/H < 1.001 && H_guess/H > 0.999)
            break
        end
    end
end