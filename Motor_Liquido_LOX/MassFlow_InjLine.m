% Esta funcion nos entregara el flujo de masa de los inyectores
% en funcion de: 1) presion del tanque de oxidante, 2) Presion del tanque de combustible, 
% 3) presion en la camara de combustion, 4) Temperatura del Oxidante,
% 5) Parametros de plometria
% 
% Para esto se calcula 
% 1)La caida de presion en la tuberia como flujos bisfasico
% 2)La densdad de flujo de masa bifasica

%[mOx mFu]
%[G_Dyer G_Omega, G_SPI, G_Omar]

function [mOx mFu]=MassFlow_InjLine(PoOx, PoFu,Pc, TOx, Pumbling)

    ConstantesFluidos

    % Parametros de la Plomeria
    DLineOx=Pumbling(1);
    DLineFu=Pumbling(2);
    LLineOx=Pumbling(3);
    LLineFu=Pumbling(4);
    Mul_Ox=Pumbling(5);
    Mul_Fu=Pumbling(6);
    kLineFu=Pumbling(7);
    kLineOx=Pumbling(8);
    Cpl=Pumbling(9);
    RolFu=Pumbling(10);
    Cd=Pumbling(11);
    AInjOx=Pumbling(12);
    AInjFu=Pumbling(13);
    ALineOx=3.14*.25*DLineOx^2;
    ALineFu=3.14*.25*DLineFu^2;

    Te1=TOx;                                                        %Temperatura de el oxidante antes del inyector
    PNos0=b1_PvapT*TOx.^3 + b2_PvapT*TOx.^2 + b3_PvapT*TOx + b4_PvapT;   % Presion de Vapor en el tanque % Sustancia
    CPR=0.7;
    Max_Omar=0;
    if Pc<PoOx && Pc<PoFu        
        %Estado antes de los inyectores
        M_Omega_Mod_guess=1;
        Xh1=0;
        Xh2=0;
        dPOx=(PoOx-Pc)*.01;
        dPFu=(PoFu-Pc)*.01;
        mFu0=0.5;
        mOx=0.5;
        for m=1:10
                                
            P1Ox=PoOx-dPOx;                                                     % Presion total antes de inyectores
            P1Fu=PoFu-dPFu;                                                     % Presion total antes de inyectores
    
            % sub enfriado supercargado
            % ¿Presion antes del inyector superior a Presion de vapor?
            if P1Ox>=PNos0
                PNos1=b1_PvapT*Te1.^3 + b2_PvapT*Te1.^2 + b3_PvapT*Te1 + b4_PvapT;   % Presion de Vapor en el tanque % Sustancia
            % Saturado
            else
               PNos1=P1Ox;                                                                  % Es establece que la presion 
               Te1 = b1_TPvap*exp(b2_TPvap*PNos1) +  b3_TPvap*exp(b4_TPvap*PNos1);      % Temperatura del oxidante depende de la temperatura, % Sustancia
            end
    
    
            % Propiedades Termodinamicas del Oxidante en el tanque
            Hl0 = (b1_HlP)*PNos1.^3 + (b2_HlP)*PNos1.^2 + (b3_HlP).*PNos1 + b4_HlP;                     % Sustancia
            Vl0 = b1_VlP*exp((b2_VlP)*PNos1) +  (b3_VlP)*exp((b4_VlP)*PNos1);                              % Sustancia
            Hg0 = (b1_HgP)*PNos1^8 + (b2_HgP)*PNos1^7 + (b3_HgP)*PNos1^6 + (b4_HgP)*PNos1^5 + ...
                  (b5_HgP)*PNos1^4 + (b6_HgP)*PNos1^3 + (b7_HgP)*PNos1^2 + (b8_HgP)*PNos1 + (b9_HgP); % Sustancia
            Sl0 = ( b1_SlP)*PNos1^3 + (b2_SlP)*PNos1^2 + (b3_SlP)*PNos1 + (b4_SlP);                            % Sustancia
            Vg0 =  b1_VgP*exp((b2_VgP )*PNos1) +  (b3_VgP)*exp((b4_VgP)*PNos1);                               % Sustancia
            Sg0 = (b1_SgP)*PNos1^7 + (b2_SgP)*PNos1^6 + (b3_SgP)*PNos1^5 + (b4_SgP)*PNos1^4 + ...
                      (b5_SgP)*PNos1^3 + (b6_SgP)*PNos1^2 + (b7_SgP )*PNos1 + (b8_SgP);                         % Sustancia
            Cpg0 = b1_Cpg*exp(b2_Cpg*TOx) +  (b3_Cpg)*exp(b4_Cpg *TOx);                                              % Sustancia
            Cpl0 = b1_Cpl*exp( b2_Cpl *TOx) + (b3_Cpl)*exp(b4_Cpl  *TOx);                                            % Sustancia
    
            
            %Calculo de la caida de Presion en tuberia y conectores
            ReOx=4*M_Omega_Mod_guess/(3.14*DLineOx*Mul_Ox);              
            AOx=( (kLineOx^1.1098)/2.8257)+(7.149/ReOx)^0.8981;
            FdOx=real(0.25*log(kLineOx/3.706-5.04*log(AOx)/ReOx)^-2);
            dPfOx=((M_Omega_Mod_guess/ALineOx)^2)*(2*FdOx*(Vl0+Xh1*(Vg0-Vl0)/2)/DLineOx + Xh1*(Vg0-Vl0));
    
            dPrOx=(dPfOx)*LLineOx;
            dPOx=dPrOx;
        
            ReFu=4*mFu0/(3.14*DLineFu*Mul_Fu);              
            AFu=( (kLineFu^1.1098)/2.8257)+(7.149/ReFu)^0.8981;
            FdFu=real(0.25*log(kLineFu/3.706-5.04*log(AFu)/ReFu)^-2);
            dPfFu=((mFu0/ALineFu)^2)*(2*FdFu*Vl0/DLineFu);
    
            dPFu=(dPfFu)*LLineFu;
            mFu0=subplus(real(Cd*AInjFu*(2*(P1Fu-Pc)*RolFu)^0.5));
            % Propiedades Termodinamicas del Oxidante antes de los Inyectores
    
            Hl1 = (b1_HlP)*PNos1.^3 + (b2_HlP)*PNos1.^2 + (b3_HlP).*PNos1 + b4_HlP;
            Vl1 = b1_VlP*exp((b2_VlP)*PNos1) +  (b3_VlP)*exp((b4_VlP)*PNos1);
            Hg1 = (b1_HgP)*PNos1^8 + (b2_HgP)*PNos1^7 + (b3_HgP)*PNos1^6 + (b4_HgP)*PNos1^5 + ...
                  (b5_HgP)*PNos1^4 + (b6_HgP)*PNos1^3 + (b7_HgP)*PNos1^2 + (b8_HgP)*PNos1 + (b9_HgP);
            Vg1 =  b1_VgP*exp((b2_VgP )*PNos1) +  (b3_VgP)*exp((b4_VgP)*PNos1);
            Cpg = b1_Cpg*exp(b2_Cpg*Te1) +  (b3_Cpg)*exp(b4_Cpg *Te1);
            Cpl = b1_Cpl*exp( b2_Cpl *Te1) + (b3_Cpl)*exp(b4_Cpl  *Te1);
    
            % Propiedades Termodinamicas del Oxidante despues de los Inyectores
    
            Vl2 = Vl1;
            Hl2 = Hl1;
            Hg2 = Hg1;  
            Vg2 = Vg1;
            Te2 = Te1;
    
            %Solo habra un cambio de fase si la presion de la camara es
            %menor que la presion de vapor. La caida de presion que primero
            %ocurrira sera la de sobrecarga
    
            if Pc<=PNos1
                Vl2 = b1_VlP*exp((b2_VlP)*Pc) +  (b3_VlP)*exp((b4_VlP)*Pc);
                Hl2 = (b1_HlP)*Pc.^3 + (b2_HlP)*Pc.^2 + (b3_HlP).*Pc + b4_HlP;
                Hg2 = (b1_HgP)*Pc^8 + (b2_HgP)*Pc^7 + (b3_HgP)*Pc^6 + (b4_HgP)*Pc^5 + ...
                      (b5_HgP)*Pc^4 + (b6_HgP)*Pc^3 + (b7_HgP)*Pc^2 + (b8_HgP)*Pc + (b9_HgP);  
                Sg2 = (b1_SgP)*Pc^7 + (b2_SgP)*Pc^6 + (b3_SgP)*Pc^5 + (b4_SgP)*Pc^4 + ...
                      (b5_SgP)*Pc^3 + (b6_SgP)*Pc^2 + (b7_SgP )*Pc + (b8_SgP);
                Sl2 = ( b1_SlP)*Pc^3 + (b2_SlP)*Pc^2 + (b3_SlP)*Pc + (b4_SlP);
                Vg2 =  b1_VgP*exp((b2_VgP )*Pc) +  (b3_VgP)*exp((b4_VgP)*Pc);
                Te2 = b1_TPvap*exp((b2_TPvap)*Pc) +  (b3_TPvap)*exp((b4_TPvap)*Pc);
            end
            
    
            Xh1 = (Hl0-Hl1)./(Hg1-Hl1);          
            Xh2 = (Hl1-Hl2)./(Hg2-Hl2);   
% ---------------------------------------------------------------------------------
            % Omega
            w=Cpl*Te1*PNos1*(1/Vl1)*((Vg1-Vl1)/(Hg1-Hl1))^2;
            CPRs=2*w/(1+2*w);
            PR=Pc/PoOx;
            PRs=PNos1/PoOx;
            CPR=CPR;
            for n=1:10
                CPR1=CPR;
                CPR=exp(-(CPR^2 + (w^2-2*w)*(1-CPR)^2 + (2*w^2)*(1-CPR))/(2*w^2));
                CPR=CPR*0.5+CPR1*0.5;
            end
            
            % Low subcooled
            if PNos1 >= CPRs*PoOx
                G_Low=subplus(real((( 2*(1-PRs) + 2*(w*PRs*log(PRs/PR)-(w-1)*(PRs-PR)) )^.5)*((PoOx/Vl1)^.5)/(w*(PRs/PR-1) + 1)));
                G_High=subplus(real(CPR*(PoOx/(w*Vl1))^.5));
                G_Omega=(PNos1/PoOx)*G_High +(1-PNos1/PoOx)*G_Low;
            
            % High subcooled 
            else
                if Pc>=PNos1
                    G_Omega=subplus(real((2*(PoOx-Pc)/Vl1)^.5));
                else
                    G_Omega=subplus(real((2*(PoOx-PNos1)/Vl1)^.5));
                end
            end         

% ---------------------------------------------------------------------------------------
            % Modelo Omar, isoentalpico considerando a entalpia de
            % fluido comprimido
            dHh=Hl1-(1-Xh2).*Hl2-Xh2.*Hg2 + (P1Ox- Pc)*((1-Xh2).*Vl2 + Xh2.*Vg2); 
            Flujo_Omar= subplus(real((2*dHh).^0.5./((1-Xh2).*Vl2 + Xh2.*Vg2)));
            
            % Si el Pressure ratio es superior a Critico
            
%             if Pc/P1Ox>=CPR
%                 G_Omar=Flujo_Omar;
%             % Si es menor se debe buscar el flujo maximo para su respectiva
%             % Po
%             else 
                DP=-0.05*10^6;
                for PcI=PoOx:DP:Pc
                    Vl2I = b1_VlP*exp((b2_VlP)*PcI) +  (b3_VlP)*exp((b4_VlP)*PcI);
                    Hl2I = (b1_HlP)*PcI.^3 + (b2_HlP)*PcI.^2 + (b3_HlP).*PcI + b4_HlP;
                    Hg2I = (b1_HgP)*PcI^8 + (b2_HgP)*PcI^7 + (b3_HgP)*PcI^6 + (b4_HgP)*PcI^5 + ...
                      (b5_HgP)*PcI^4 + (b6_HgP)*PcI^3 + (b7_HgP)*PcI^2 + (b8_HgP)*PcI + (b9_HgP);  
                    Vg2I =  b1_VgP*exp((b2_VgP )*PcI) +  (b3_VgP)*exp((b4_VgP)*PcI);
                    Xh2I =       (Hl1-Hl2I)./(Hg2I-Hl2I);   
 
                    dHhI=subplus(Hl1-(1-Xh2I).*Hl2I-Xh2I.*Hg2I + (P1Ox- PcI)*((1-Xh2I).*Vl2I + Xh2I.*Vg2I));
                    Flujo_OmarI= (2*dHhI).^0.5./((1-Xh2I).*Vl2I + Xh2I.*Vg2I);
                    if Flujo_OmarI>=Max_Omar
                        Max_Omar=Flujo_OmarI;
                        G_Omar=Flujo_OmarI;
                    end
  
                end
%             end
% --------------------------------------------------------------------------------    
            % Dyer
            G_SPI=subplus(2*(PoOx-Pc)/Vl1).^.5;       
            k_Dyer=subplus((P1Ox-Pc)/(PNos1-Pc))^.5;
            if Pc>=PNos1
                G_Dyer=G_SPI;
            else
                G_Dyer=k_Dyer*G_SPI/(k_Dyer+1)+G_Omar/(k_Dyer+1);
            end  
% -------------------------------------------------------------------------------
            % Omega Modificado
            t=3.1416*.5*(P1Ox-Pc)/(P1Ox-P1Ox*CPR);
            t1=sin(t)^.5;
            t2=1-sin(t)^.5;
            if Pc>=P1Ox*CPR            
                M_Omega_Mod=G_Omega*t1+ G_Dyer*t2;
            else
                M_Omega_Mod=G_Omega;
            end
            M_Omega_Mod=Cd*AInjOx*M_Omega_Mod;
            M_Omega_Mod_guess=real(M_Omega_Mod*0.5 + M_Omega_Mod_guess*0.5);
        end
    else
        G_Omar=0;   
        G_SPI=0;
        G_Dyer=0;
        G_Omega=0;
        M_Omega_Mod=0;
        mFu0=0;
    end
    mFu=mFu0;
    mOx=G_Omar*Cd*AInjOx;
end