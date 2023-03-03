clear all
close all

ConstantesFluidos



%Parametros de la Tuberia
LLineOx=1;
WLineOx=0.889*10^-3;
DLineOx=0.0254*4/8-2*WLineOx;
ALineOx=3.14*.25*DLineOx^2;
K_Cooper=50*WLineOx
kLineOx=3*10^-6/DLineOx; % Rugosidad relativa de la tuberia de cobre o acero inoxidable flexible "Cooper tubing"

% Parametros de lso inyectores
LInjOx=10*10^-3
DInjOx=2*10^-3
NInjOx=20;
AInjOx=NInjOx*3.14*.25*DInjOx^2

T_Atm=290;
hPipeAat=50;


Mul_Ox=  0.000077;

%Estado del NOS en el Tanque
Psup=0.0*10^6;    % Presion de supercarga del N2 al NOS
Te0=270;        % Temperatura de Saturacion del NOS
Te1=Te0;
PNos0=2.792*Te0^3 + -1527*Te0^2 + 2.898e+05*Te0 + -1.899e+07 % Presion de Vapor en el tanque
CPR=.7

% Limitar el rango del analisis, ya que este programa es solo para analisar
% los modelos de flujo bifasico
Pc_min=1*10^6
d_P=-0.1*10^5;
Po= PNos0 + Psup

% Inicializacion del tamaño de lso arreglos
LengthData=size(Po:d_P:Pc_min)
G_Omega=zeros([1,LengthData]);
G_Omar=G_Omega;
G_SPI=G_Omega;
G_Dyer=G_Omega;
G_High=G_Omega;
G_Low=G_Omega;
Xh1=G_Omega;
Xh2=G_Omega;
t=zeros([LengthData,1]);
Te1=Te0*ones([1,LengthData])
Hl0=t;Hg0=t;
Vl0=t;Vg0=t;
Sl0=t;Sg0=t;

Hl1=t;Hg1=t;
Vl1=t;Vg1=t;
Sl1=t;Sg1=t;

Hl2=t;Hg2=t;
Vl2=t;Vg2=t;
Sl2=t;Sg2=t;

t1=t;
t2=t;
dP(1)=0;

M_Omega_Mod_guess=1;
Max_Omar=0;
i=1;
j=1;

%%
    for Pc=Po:d_P:Pc_min  
        if Pc<Po         
            %Estado antes de los inyectores
            M_Omega_Mod_guess=1;
            dP(j)=(Po-Pc)*.01;
            dPr(j)=dP(j);
            for m=1:10
                                
                P1(i,j)=Po-dP(j);
                
                % sub enfriado supercargado
                if P1(i,j)>=PNos0
                   PNos1(i,j)=2.792*Te1(i,j)^3 + -1527*Te1(i,j)^2 + 2.898e+05*Te1(i,j) + -1.899e+07; % Presion de Vapor en el tanque                
                % Saturado
                else
                   PNos1(i,j)=P1(i,j);
                   Te1(i,j) = 268.2*exp((2.067e-08)*PNos1(i,j)) +  (-64.09)*exp((-5.139e-07)*PNos1(i,j));
                end

    
                % Propiedades Termodinamicas del Oxidante en el tanque
                Hl0(j) = (5.842e-16)*PNos1(i,j).^3 + (-8.466e-09)*PNos1(i,j).^2 + (0.06278).*PNos1(i,j) + -4.315e+05;
                Vl0(j) = 0.000858*exp((8.077e-08)*PNos1(i,j)) +  (9.149e-12)*exp((2.5e-06)*PNos1(i,j));
                Hg0(j) = (-1.223e-47)*PNos1(i,j)^8 + (3.658e-40)*PNos1(i,j)^7 + (-4.626e-33)*PNos1(i,j)^6 + (3.221e-26)*PNos1(i,j)^5 + ...
                      (-1.346e-19)*PNos1(i,j)^4 + (3.452e-13)*PNos1(i,j)^3 + (-5.306e-07)*PNos1(i,j)^2 + (0.4495)*PNos1(i,j) + (-2.311e+05);
                Sl0(j) = ( 2.44e-18)*PNos1(i,j)^3 + (-3.679e-11)*PNos1(i,j)^2 + (0.0002576)*PNos1(i,j) + (-248.2);
                Vg0(j) =  0.09481*exp((-1.873e-06 )*PNos1(i,j)) +  (0.03302)*exp((-3.344e-07)*PNos1(i,j));
                Sg0(j) = (-1.296e-43)*PNos1(i,j)^7 + (3.411e-36)*PNos1(i,j)^6 + (-3.702e-29)*PNos1(i,j)^5 + (2.141e-22)*PNos1(i,j)^4 + ...
                          (-7.116e-16)*PNos1(i,j)^3 + (1.366e-09)*PNos1(i,j)^2 + (-0.001477 )*PNos1(i,j) + (1919);
                Cpg0 = 308.4*exp(0.005439*Te0) +  (6.942e-13)*exp(0.121 *Te0);
                Cpl0 = 665*exp( 0.004254 *Te0) + (8.183e-10)*exp(0.09463  *Te0);
    
                
     
                %Calculo de la caida de Presion en tuberia y conectores
                Re=4*M_Omega_Mod_guess/(3.14*DLineOx*Mul_Ox);              
                A=( (kLineOx^1.1098)/2.8257)+(7.149/Re)^0.8981;
                Fd=0.25*log(kLineOx/3.706-5.04*log(A)/Re)^-2;
                dPf=((M_Omega_Mod_guess/ALineOx)^2)*(2*Fd*(Vl0(j)+Xh1(i,j)*(Vg0(j)-Vl0(j))/2)/DLineOx + Xh1(i,j)*(Vg0(j)-Vl0(j)));

                dPr(j)=(dPf)*LLineOx;
                dP(j)=dPr(j);

                Cd=0.66;
                % Propiedades Termodinamicas del Oxidante antes de los Inyectores
    
                Hl1(j) = (5.842e-16)*PNos1(i,j).^3 + (-8.466e-09)*PNos1(i,j).^2 + (0.06278).*PNos1(i,j) + -4.315e+05;
                Vl1(j) = 0.000858*exp((8.077e-08)*PNos1(i,j)) +  (9.149e-12)*exp((2.5e-06)*PNos1(i,j));
                Hg1(j) = (-1.223e-47)*PNos1(i,j)^8 + (3.658e-40)*PNos1(i,j)^7 + (-4.626e-33)*PNos1(i,j)^6 + (3.221e-26)*PNos1(i,j)^5 + ...
                      (-1.346e-19)*PNos1(i,j)^4 + (3.452e-13)*PNos1(i,j)^3 + (-5.306e-07)*PNos1(i,j)^2 + (0.4495)*PNos1(i,j) + (-2.311e+05);
                Sl1(j) = ( 2.44e-18)*PNos1(i,j)^3 + (-3.679e-11)*PNos1(i,j)^2 + (0.0002576)*PNos1(i,j) + (-248.2);
                Vg1(j) =  0.09481*exp((-1.873e-06 )*PNos1(i,j)) +  (0.03302)*exp((-3.344e-07)*PNos1(i,j));
                Sg1(j) = (-1.296e-43)*PNos1(i,j)^7 + (3.411e-36)*PNos1(i,j)^6 + (-3.702e-29)*PNos1(i,j)^5 + (2.141e-22)*PNos1(i,j)^4 + ...
                          (-7.116e-16)*PNos1(i,j)^3 + (1.366e-09)*PNos1(i,j)^2 + (-0.001477 )*PNos1(i,j) + (1919);
                Cpg = 308.4*exp(0.005439*Te1(i,j)) +  (6.942e-13)*exp(0.121 *Te1(i,j));
                Cpl = 665*exp( 0.004254 *Te1(i,j)) + (8.183e-10)*exp(0.09463  *Te1(i,j));
    
                % Propiedades Termodinamicas del Oxidante despues de los Inyectores
    
                Vl2(j) = Vl1(j);
                Hl2(j) = Hl1(j);
                Hg2(j) = Hg1(j);  
                Sg2(j) = Sg1(j);
                Sl2(j) = Sl1(j);
                Vg2(j) = Vg1(j);
                Te2(i,j) = Te1(i,j);
                dSl2_dP = 3*( 2.44e-18)*PNos1(i,j)^2 + 2*(-3.679e-11)*PNos1(i,j) + (0.0002576);
    
                %Solo habra un cambio de fase si la presion de la camara es
                %menor que la presion de vapor. La caida de presion que primero
                %ocurrira sera la de sobrecarga

                if Pc<=PNos1(i,j)
                    Vl2(j) = 0.000858*exp((8.077e-08)*Pc) +  (9.149e-12)*exp((2.5e-06)*Pc);
                    Hl2(j) = (5.842e-16)*Pc.^3 + (-8.466e-09)*Pc.^2 + (0.06278).*Pc + -4.315e+05;
                    Hg2(j) = (-1.223e-47)*Pc^8 + (3.658e-40)*Pc^7 + (-4.626e-33)*Pc^6 + (3.221e-26)*Pc^5 + ...
                          (-1.346e-19)*Pc^4 + (3.452e-13)*Pc^3 + (-5.306e-07)*Pc^2 + (0.4495)*Pc + (-2.311e+05);  
                    Sg2(j) = (-1.296e-43)*Pc^7 + (3.411e-36)*Pc^6 + (-3.702e-29)*Pc^5 + (2.141e-22)*Pc^4 + ...
                          (-7.116e-16)*Pc^3 + (1.366e-09)*Pc^2 + (-0.001477 )*Pc + (1919);
                    Sl2(j) = ( 2.44e-18)*Pc^3 + (-3.679e-11)*Pc^2 + (0.0002576)*Pc + (-248.2);
                    Vg2(j) =  0.09481*exp((-1.873e-06 )*Pc) +  (0.03302)*exp((-3.344e-07)*Pc);
                    Te2(i,j) = 268.2*exp((2.067e-08)*Pc) +  (-64.09)*exp((-5.139e-07)*Pc);
                    dSl2_dP = 3*( 2.44e-18)*Pc^2 + 2*(-3.679e-11)*Pc + (0.0002576);
                end
                

                Xh1(i,j) = (Hl0(j)-Hl1(j))./(Hg1(j)-Hl1(j));          
                Xh2(i,j) =       (Hl1(j)-Hl2(j))./(Hg2(j)-Hl2(j));   
                % Omega
                w=Cpl*Te1(i,j)*PNos1(i,j)*(1/Vl1(j))*((Vg1(j)-Vl1(j))/(Hg1(j)-Hl1(j)))^2;
                CPRs=2*w/(1+2*w);
                PR=Pc/Po;
                PRs=PNos1(i,j)/Po;
                CPR=CPR;
                for n=1:10
                    CPR1=CPR;
                    CPR=exp(-(CPR^2 + (w^2-2*w)*(1-CPR)^2 + (2*w^2)*(1-CPR))/(2*w^2));
                    CPR=CPR*0.5+CPR1*0.5;
                end
                
                % Low subcooled
                if PNos1(i,j) >= CPRs*Po
                    G_Low(i, j)=(( 2*(1-PRs) + 2*(w*PRs*log(PRs/PR)-(w-1)*(PRs-PR)) )^.5)*((Po/Vl1(j))^.5)/(w*(PRs/PR-1) + 1);
                    G_High(i,j)=CPR*(Po/(w*Vl1(j)))^.5;
                    G_Omega(i,j)=(PNos1(i,j)/Po)*G_High(i,j) +(1-PNos1(i,j)/Po)*G_Low(i,j);
                
                % High subcooled 
                else
                    if Pc>=PNos1(i,j)
                    G_Omega(i,j)=(2*(Po-Pc)/Vl1(j))^.5;
                    else
                    G_Omega(i,j)=(2*(Po-PNos1(i,j))/Vl1(j))^.5;
                    end
                end         
                
                % Dyer
                [m P1(i,j) Pc];
                G_SPI(i,j)=(2*(Po-Pc)/Vl1(j)).^.5;       
                k_Dyer=((P1(i,j)-Pc)/(PNos1(i,j)-Pc))^.5;
                if Pc>=PNos1(i,j)
                    G_Dyer(i,j)=G_SPI(i,j);
                else
                    G_Dyer(i,j)=k_Dyer*G_SPI(i,j)/(k_Dyer+1)+G_Omar(i,j)/(k_Dyer+1);
                end
    
                % Modelo Omar, isoentalpico considerando a entalpia de
                % fluido comprimido
                dHh=Hl1(j)-(1-Xh2(i,j)).*Hl2(j)-Xh2(i,j).*Hg2(j) + (P1(i,j)- Pc)*((1-Xh2(i,j)).*Vl2(j) + Xh2(i,j).*Vg2(j)) 
                Flujo_Omar= (2*dHh).^0.5./((1-Xh2(i,j)).*Vl2(j) + Xh2(i,j).*Vg2(j));
                if Flujo_Omar>Max_Omar
                    Max_Omar=Flujo_Omar;
                end
                G_Omar(i,j)=Max_Omar;

                % Omega Modificado
                t(j)=3.1416*.5*(P1(i,j)-Pc)/(P1(i,j)-P1(i,j)*CPR);
                t1(j)=sin(t(j))^.5;
                t2(j)=1-sin(t(j))^.5;
                if Pc>=P1(i,j)*CPR            
                    M_Omega_Mod(i,j)=G_Omega(i,j)*t1(j)+ G_Dyer(i,j)*t2(j);
                else
                    M_Omega_Mod(i,j)=G_Omega(i,j);
                end
                M_Omega_Mod(i,j)=Cd*AInjOx*M_Omega_Mod(i,j);
                M_Omega_Mod_guess=real(M_Omega_Mod(i,j)*0.5 + M_Omega_Mod_guess*0.5);
            end
        else
            G_Omar(i,j)=0;   
            G_SPI(i,j)=0;
            G_Dyer(i,j)=0;
        end
        j=j+1;
   
    end


Pc=Po:d_P:Pc_min;
[Po,Pc]=meshgrid(Po,Pc);

figure (1)
plot( Pc, Cd*AInjOx*M_Omega_Mod(1,:), Pc, Cd*AInjOx*G_Omega(1,:), Pc, Cd*AInjOx*G_Omar(1,:), Pc, Cd*AInjOx*G_SPI(1,:), Pc, Cd*AInjOx*G_Dyer(1,:), Pc, Cd*AInjOx*G_High(1,:), Pc, Cd*AInjOx*G_Low(1,:))
legend('OmegaMod','Omega', 'Omar', 'SPI', 'Dyer', 'High', 'Low')
title("Mass Flow")

figure (2)
plot(Pc, Te2(1,:), Pc, Te1(1,:), Pc, Te0*ones([size(Te1(1,:)),1]))
title('Temperatura NOS Camara')
legend('T2', 'T1', 'T0')

figure (3)
plot(Pc, P1(1,:), Pc, PNos1(1,:), Pc, dPr, Pc, Pc, Pc, Po, Pc, P1(1,:)+dPr, Pc, dP)
title('Presiones')
legend('PInj', 'PNos1', 'dPr', 'Pc', 'Po' ,'PInj+dPr', 'dP')

figure (4)
plot(Pc, Xh1(1,:), Pc, Xh2(1,:))
title('Calidad Equilibrio')
legend('XPipe', 'XChamber')


figure (5)
plot(Pc, Hl0, Pc, Hg0, Pc, Hl1, Pc, Hg1, Pc, Hl2, Pc, Hg2)
title('Entalpias')
legend ('l0', 'g0', 'l1', 'g1', 'l2', 'g2')

figure (6)
plot(Pc, Vl0, Pc, Vg0, Pc, Vl1, Pc, Vg1, Pc, Vl2, Pc, Vg2)
title('Volumen espec')
legend ('l0', 'g0', 'l1', 'g1', 'l2', 'g2')
