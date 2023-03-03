syms P1 P2 Hl1 Hl2 Hg1 Hg2 Sl1 Sl2 Sg1 Sg2 Vg1 Vg2 Vl1 Vl2 Te1 Te2 X1 X2 MHEM
eq1 = Hl1 == -3.826e-06*P1^2 +  0.0516*P1    + -423.4
eq11= Hl2 == -3.826e-06*P2^2 +  0.0516*P2    + -423.4

eq2 = Hg1 == -1.869e-06*P1^2 +  0.007673*P1  + -77.25
eq21= Hg2 == -1.869e-06*P2^2 +  0.007673*P2  + -77.25

eq31= Sg2 ==  5.786e-09*P2^2 + -9.386e-05*P2 +  1.344

eq4 = Sl1 == -1.736e-08*P1^2 +  0.0002114*P1 + -0.2151
eq41= Sl2 == -1.736e-08*P2^2 +  0.0002114*P2 + -0.2151

eq5 = Vg1 ==  2.333e-09*P1^2 + -2.016e-05*P1 + 0.05129
eq51= Vg2 ==  2.333e-09*P2^2 + -2.016e-05*P2 + 0.05129

eq6 = Vl1 ==                   7.758e-08*P1 +  0.0008631
eq61= Vl2 == 0.00108      %            7.758e-08*P2 +  0.0008631

eq8=  Te1 == -1.923e-06*P1^2 +  0.02518*P1   + 213.5
eq81= Te2 == -1.923e-06*P2^2 +  0.02518*P2   + 213.5

eq9= X2  == (Sl1-Sl2)/(Sg2-Sl2)
eq10= MHEM== ((2*(Hl1-(1-X2)*Hl2-X2*Hg2))^.5)/( (1-X2)*Vl2 + X2*Vg2)

EQ=[eq1, eq11, eq21, eq31, eq4, eq41, eq51, eq61, eq9, eq10]
X =[Hl1,  Hl2,  Hg2,  Sg2, Sl1,  Sl2,  Vg2,  Vl2,  X2, MHEM]
Sol=solve (EQ, X)

vpa(Sol.MHEM)
vpa(Sol.X2)

%%
clear all
close all

% M_HEM1=-(1.0.*(0.1032.*P11 - 2.0.*((- 3279211285704681875.0.*P11.^2 + 39932330979145717760000.0.*P11 + 3279211285704681875.0.*P21.^2 - 39932330979145717760000.0.*P21)./(4372155784500032500.0.*P21.^2 - 57661983702431513600000.0.*P21 + 294505663337682554106937344.0) - 1.0).*(0.0000038260000000000002990871439401133.*P21.^2 - 0.0516.*P21 + 423.4) - 0.0000076520000000000005981742878802265.*P11.^2 + (2.0.*(0.0000018690000000000000202445482253011.*P21.^2 - 0.0076730000000000001328381848964.*P21 + 77.25).*(- 3279211285704681875.0.*P11.^2 + ...
%     39932330979145717760000.0.*P11 + 3279211285704681875.0.*P21.^2 - 39932330979145717760000.0.*P21))./(4372155784500032500.0.*P21.^2 - 57661983702431513600000.0.*P21 + 294505663337682554106937344.0) - 846.8).^(1./2))./((0.000000077579999999999997209336520618411.*P21 + 0.00086310000000000000001249000902703).*((- 3279211285704681875.0.*P11.^2 + 39932330979145717760000.0.*P11 + 3279211285704681875.0.*P21.^2 - 39932330979145717760000.0.*P21)./(4372155784500032500.0.*P21.^2 - 57661983702431513600000.0.*P21 + 294505663337682554106937344.0) - 1.0) ...
%     - (1.0.*(0.000000002333000000000000112008933215741.*P21.^2 - 0.000020160000000000000104161471115027.*P21 + 0.05129).*(- 3279211285704681875.0.*P11.^2 + 39932330979145717760000.0.*P11 + 3279211285704681875.0.*P21.^2 - 39932330979145717760000.0.*P21))./(4372155784500032500.0.*P21.^2 - 57661983702431513600000.0.*P21 + 294505663337682554106937344.0))
% X1=(- 3279211285704681875.0.*P11.^2 + 39932330979145717760000.0.*P11 + 3279211285704681875.0.*P21.^2 - 39932330979145717760000.0.*P21)./(4372155784500032500.0.*P21.^2 - 57661983702431513600000.0.*P21 + 294505663337682554106937344.0)

M_HEM2=zeros([1,9]);
M_ERM=M_HEM2;
M_Omar=M_HEM2;
M_SPI=M_HEM2;
M_SPIg=M_HEM2;
M_Dyer=M_HEM2;
M_Omega=M_HEM2;

M_RLM=M_HEM2;
M_HFSM=M_HEM2;

Xs=M_HEM2;
Xh=M_HEM2;

Vl2m=M_HEM2;
Vmix=M_HEM2;
i=1;
L=0.010;
D=0.001;
k=2*10^-6/D;
M_Max=0;
gamma=1.31
R=3.14/0.044;

Psup=0.5*10^6;
Te1=270;
PNos0=2.792*Te1^3 + -1527*Te1^2 + 2.898e+05*Te1 + -1.899e+07

Pc_min=1*10^6
Po= PNos0 + Psup
CPR=.7

for Po=Po
    j=1;
    M_Max=0;

    for Pc=Po:-0.1*10^5:Pc_min  
        if Pc<Po 

            dP=.1*10^6

            P1=Po-dP

            if P1<=PNos0
               PNos1=P1
            else
               PNos1=PNos0
            end

            Hl1 = (5.842e-16)*PNos1.^3 + (-8.466e-09)*PNos1.^2 + (0.06278).*PNos1 + -4.315e+05;
            Vl1 = 0.000858*exp((8.077e-08)*PNos1) +  (9.149e-12)*exp((2.5e-06)*PNos1);
            Hg1 = (-1.223e-47)*PNos1^8 + (3.658e-40)*PNos1^7 + (-4.626e-33)*PNos1^6 + (3.221e-26)*PNos1^5 + ...
                  (-1.346e-19)*PNos1^4 + (3.452e-13)*PNos1^3 + (-5.306e-07)*PNos1^2 + (0.4495)*PNos1 + (-2.311e+05);
            Sl1 = ( 2.44e-18)*PNos1^3 + (-3.679e-11)*PNos1^2 + (0.0002576)*PNos1 + (-248.2);
            Vg1 =  0.09481*exp((-1.873e-06 )*PNos1) +  (0.03302)*exp((-3.344e-07)*PNos1);
            Sg1 = (-1.296e-43)*PNos1^7 + (3.411e-36)*PNos1^6 + (-3.702e-29)*PNos1^5 + (2.141e-22)*PNos1^4 + ...
                      (-7.116e-16)*PNos1^3 + (1.366e-09)*PNos1^2 + (-0.001477 )*PNos1 + (1919);
            Cpg = 308.4*exp(0.005439*Te1) +  (6.942e-13)*exp(0.121 *Te1);
            Cpl = 665*exp( 0.004254 *Te1) + (8.183e-10)*exp(0.09463  *Te1);

            Vl2 = Vl1;
            Hl2 = Hl1;
            Hg2 = Hg1;  
            Sg2 = Sg1;
            Sl2 = Sl1;
            Vg2 = Vg1;
            Te2(i,j) = Te1;
            dSl2_dP = 3*( 2.44e-18)*PNos1^2 + 2*(-3.679e-11)*PNos1 + (0.0002576);

            %Solo habra un cambio de fase si la presion de la camara es
            %menor que la presion de vapor. La caida de presion que primero
            %ocurrira sera la de sobrecarga

            if Pc<=PNos1
                Vl2 = 0.000858*exp((8.077e-08)*Pc) +  (9.149e-12)*exp((2.5e-06)*Pc);
                Hl2 = (5.842e-16)*Pc.^3 + (-8.466e-09)*Pc.^2 + (0.06278).*Pc + -4.315e+05;
                Hg2 = (-1.223e-47)*Pc^8 + (3.658e-40)*Pc^7 + (-4.626e-33)*Pc^6 + (3.221e-26)*Pc^5 + ...
                      (-1.346e-19)*Pc^4 + (3.452e-13)*Pc^3 + (-5.306e-07)*Pc^2 + (0.4495)*Pc + (-2.311e+05);  
                Sg2 = (-1.296e-43)*Pc^7 + (3.411e-36)*Pc^6 + (-3.702e-29)*Pc^5 + (2.141e-22)*Pc^4 + ...
                      (-7.116e-16)*Pc^3 + (1.366e-09)*Pc^2 + (-0.001477 )*Pc + (1919);
                Sl2 = ( 2.44e-18)*Pc^3 + (-3.679e-11)*Pc^2 + (0.0002576)*Pc + (-248.2);
                Vg2 =  0.09481*exp((-1.873e-06 )*Pc) +  (0.03302)*exp((-3.344e-07)*Pc);
                Te2(i,j) = 268.2*exp((2.067e-08)*Pc) +  (-64.09)*exp((-5.139e-07)*Pc);
                dSl2_dP = 3*( 2.44e-18)*Pc^2 + 2*(-3.679e-11)*Pc + (0.0002576);
            end
            
            % Omega
            w=Cpl*Te1*PNos1*(1/Vl1)*((Vg1-Vl1)/(Hg1-Hl1))^2
            CPRs=2*w/(1+2*w)
            PR=Pc/Po
            PRs=PNos1/Po
            CPR=CPR
            for n=1:10
                CPR1=CPR;
                CPR=exp(-(CPR^2 + (w^2-2*w)*(1-CPR)^2 + (2*w^2)*(1-CPR))/(2*w^2));
                CPR=CPR*0.5+CPR1*0.5;
                %(-p2-(p2^2-4*p1*p3)^.5)/(2*p1)
            end
            % Low subcooled
            if PNos1 >= CPRs*Po
            
            G_Low=(( 2*(1-PRs) + 2*(w*PRs*log(PRs/PR)-(w-1)*(PRs-PR)) )^.5)*((Po/Vl1)^.5)/(w*(PRs/PR-1) + 1)
            G_High=CPR*(Po/(w*Vl1))^.5
            M_Omega(i,j)=(PNos1/Po)*G_High +(1-PNos1/Po)*G_Low;
            % High subcooled 
            else
                if Pc>=PNos1
                M_Omega(i,j)=(2*(Po-Pc)/Vl1)^.5
                else
                M_Omega(i,j)=(2*(Po-PNos1)/Vl1)^.5
                end
            end

            Mu_l=  0.000077;
            Xs(i,j)  = (Sl1-Sl2)./(Sg2-Sl2);
            Xh(i,j) = (Hl1-Hl2)./(Hg2-Hl2);          
            
            dHs=Hl1-(1-Xs(i,j)).*Hl2-Xs(i,j).*Hg2 + (Po-Pc)*((1-Xs(i,j)).*Vl2 + Xs(i,j).*Vg2);
            dHh=Hl1-(1-Xh(i,j)).*Hl2-Xh(i,j).*Hg2 + (Po-Pc)*((1-Xh(i,j)).*Vl2 + Xh(i,j).*Vg2);

            M_ERM(i,j)=(Hg2-Hl2)/((Vg2-Vl2)*(Te2(i,j)*Cpl)^(.5));
           
            M_RLM_guess1= 1000;
            M_RLM_guess2= 1000;
            for n=1:10
                M_RLM_guess1= M_RLM_guess1*.5 + M_RLM_guess2*.5;
                Re=M_RLM_guess1*D/Mu_l;               
                A=( (k^1.1098)/2.8257)+(7.149/Re)^0.8981;
                fd=0.25*log(k/3.706-5.04*log(A)/Re)^-2;
                K=fd*L/D;
                N=(Vl1/(2*(Po-Pc)*Cpl*Te1*K^2) )*((Hg2-Hl2)/(Vg2-Vl2))^2 + L/.1;
    
                M_RLM_guess2=(Hg2-Hl2)/((Vg2-Vl2)*(N*Te2(i,j)*Cpl)^(.5));;
                Error=(M_RLM_guess1-M_RLM_guess2)/M_RLM_guess1;

            end 
            
            M_RLM(i,j)=M_RLM_guess2;
            M_HEM2(i,j)= (2*dHs).^0.5./((1-Xs(i,j)).*Vl2 + Xs(i,j).*Vg2);
            M_Omar(i,j)= (2*dHh).^0.5./((1-Xh(i,j)).*Vl2 + Xh(i,j).*Vg2);

            M_SPI(i,j)=(2*(Po-Pc)/Vl2).^.5;
            M_SPIg(i,j)=Po*((2*gamma/(R*Te1*(gamma-1)))^.5)*((Pc/Po)^(2/gamma)-(Pc/Po)^(1+1/gamma))^.5;
            
            k_Dyer=((Po-Pc)/(PNos1-Pc))^.5;
            if Pc>=PNos1
                M_Dyer(i,j)=M_SPI(i,j);

            else
                M_Dyer(i,j)=k_Dyer*M_SPI(i,j)/(k_Dyer+1)+M_HEM2(i,j)/(k_Dyer+1);
            end

            Ns=satlin(Xs(i,j)./0.14);
            CritPressRatio=(2/(gamma+1))^(gamma/(gamma-1));
            
            Vmix_HFM=(1-Xs(i,j))*Vl1+Xs(i,j)*Vg1*CritPressRatio^(1/gamma);
            M_HFM(i,j)=(2*Xs(i,j)*Vg1*Po*(gamma/(gamma-1))^(1-CritPressRatio^((gamma-1)/gamma)))^.5/Vmix_HFM;
            M_HFSM(i,j)=((Vg2-Vl1)*(1-Xs(i,j))*Ns*dSl2_dP/(Sg2-Sl2))^-.5;
            if M_Max<M_Omar(i,j)
                M_Max=M_Omar(i,j);
            else
                M_Omar(i,j)=M_Max;  
            end
               
        else
            M_HEM2(i,j)=0;
            Xs(i,j)
            M_Omar(i,j)=0;        

        end
        j=j+1;
    end
    i=i+1;
end


Pc=Po:-0.1*10^5:Pc_min;
[Po,Pc]=meshgrid(Po,Pc);

% figure(1)
% plot3( Po,Pc, M_Omar, Po, Pc, M_SPI)
% xlabel('Po')
% ylabel('Pc')
% title('Flujo Omar & Flujo SPI')
% 
% figure(2)
% plot3 (Po,Pc, X2)
% xlabel('Po')
% ylabel('Pc')
% title('Calidad')
% 
% figure(3)
% plot3 (Po,Pc, M_RLM, Po,Pc, M_HFSM)
% xlabel('Po')
% ylabel('Pc')
% title('Rate Length Model & Hause Fause')
% 
% figure(4)
% plot3( Po,Pc, M_HEM2, Po, Pc, M_ERM)
% xlabel('Po')
% ylabel('Pc')
% title('HEM & ERM')

figure (1)
p1=plot( Pc, M_HEM2(1,:), Pc, M_SPI(1,:), Pc, M_SPIg(1,:), Pc, M_Omar(1,:), Pc, M_Dyer(1,:), Pc, M_Omega(1,:))
p1(1).LineWidth=2;
p1(2).LineWidth=2;
p1(3).LineWidth=2;
p1(4).LineWidth=2;
p1(5).LineWidth=2;
p1(6).LineWidth=2;
legend('HEM', 'SPI-L', 'SPI-G',  'Omar', 'Dyer', 'Omega')
title("Mass Flux")

figure (2)
plot(Pc, Te2(1,:))
title('Temperatura NOS Camara')

figure (3)
plot(Pc, M_HEM2(1,:))
title('HEM')

figure (4)
plot(Pc, Xs(1,:), Pc, Xh(1,:))
title('Calidad Equilibrio')
legend('S', 'H')

