clear all
%Parametros Objetivo
Impulse_o=50*1000;
Thrust=2500;
%Parametros de Combustible
Oxidizer = 'N2O';   % o 'O2'
Fuel = 'C2H6O';     % o 'n-C8H18',
Pressurant='N2';    % o 'He'
Isp_o=220;
OxFR=4;
R=8.314;        % Constnte universal de los gases
R_Press=296;    %N2,Helio 1030
Gamma_Press=1.4;%N2, Helio 1.667
R_Ox=188;       %Nitrous Oxide 
R_Fuel=181;     %Ethanol
mol_Press=0.014;  %Nitrogen, 0.002 Helium
mol_Fuel=0.046;   %Ethanol
mol_Ox=0.044;     %NOS, Oxygen 0.016
P_atm=1*10^5;
%Parametros de Tanque
T_Ox_Tank_o=263;                %Temperatura inicial del tanque de oxidante
T_Ox_Tank_u=240;                %Temperatura Final del tanque de oxidante

P_Tank_Ox=5*10^6;               %Presion total del tanque de oxidante

T_Fuel_Tank_o=300;                %Temperatura total inicial del tanque de combustible
T_Fuel_Tank_u=260;              %Temperatura Final del tanque de oxidante
P_Fuel_Tank=5*10^6;             %Presion total del tanque de combustible

T_Press_o=300;                    %Temperatura inicial del tanque de Pressurizante
T_Press_u=260;
P_Press_o=4500*6895;                %Presion Inicial del tanque de Pressurizante
P_Press_u=10*10^6;                  %Presion final del tanque de Pressurizante

%Parametros de Diseño Hidraulico
Diam_Ox_Tank=6*0.0254;          %Diametro comercial del tanque
Diam_Fuel_Tank=6*0.0254;        %Diametro comercial del tanque

%% Definir estado del tanque como mezcla
P_Ox_Vap=py.CoolProp.CoolProp.PropsSI('P','T',T_Ox_Tank_o,'Q',0,Oxidizer);         
D_Ox_L = py.CoolProp.CoolProp.PropsSI('D','T',T_Ox_Tank_o,'Q',0,Oxidizer);
D_Ox_G_u = py.CoolProp.CoolProp.PropsSI('D','T',T_Ox_Tank_u,'Q',1,Oxidizer) ; 
D_Ox_G_o = py.CoolProp.CoolProp.PropsSI('D','T',T_Ox_Tank_o,'Q',1,Oxidizer) ; 


P_Fuel_Vap=py.CoolProp.CoolProp.PropsSI('P','T',T_Fuel_Tank_o,'Q',0,Fuel);
D_Fuel_L = py.CoolProp.CoolProp.PropsSI('D','T',T_Fuel_Tank_o,'Q',0,Fuel);
D_Fuel_G = py.CoolProp.CoolProp.PropsSI('D','T',T_Fuel_Tank_o,'Q',1,Fuel);


%Ley de Roault & Ley de Dalton
P_Ox=P_Ox_Vap;              %Solo al incio, cuando las fases liquida vapor estan en equilibrio
P_Fuel=P_Fuel_Vap;          %Solo al incio, cuando las fases liquida vapor estan en equilibrio
P_Press_Ox=P_Tank_Ox-P_Ox;
P_Press_Fuel=P_Fuel_Tank-P_Fuel;

%R_Humidity=P_Ox/P_Ox_Vap;

%
D_Press_Ox_u     =   py.CoolProp.CoolProp.PropsSI('D','T',T_Ox_Tank_u,  'P',P_Press_Ox,  Pressurant)
D_Press_Ox_o     =   py.CoolProp.CoolProp.PropsSI('D','T',T_Ox_Tank_o,  'P',P_Press_Ox,  Pressurant)
D_Press_Fuel   =   py.CoolProp.CoolProp.PropsSI('D','T',T_Fuel_Tank_u,'P',P_Press_Fuel,Pressurant)
D_Press_Tank_o  =  py.CoolProp.CoolProp.PropsSI('D','T',T_Press_o,    'P',P_Press_o,   Pressurant)
D_Press_Tank_u  =  py.CoolProp.CoolProp.PropsSI('D','T',T_Press_u,    'P',P_Press_u,   Pressurant)

%Entalpias Iniciales del Oxidante y del presurizante
hmolar_Press_Press_o= py.CoolProp.CoolProp.PropsSI('Hmolar','T',T_Press_o,'P',P_Press_o, Pressurant)
h_Press_Press_o= py.CoolProp.CoolProp.PropsSI('H','T',T_Press_o,'P',P_Press_o, Pressurant)
umolar_Ox_Ox_o      = py.CoolProp.CoolProp.PropsSI('Umolar','Q',1,'P',P_Ox_Vap,  Oxidizer)
umolar_Fuel_o       = py.CoolProp.CoolProp.PropsSI('Umolar','T',T_Fuel_Tank_o ,'Q',1,  Fuel)

Mass_Prop=Impulse_o/(Isp_o*9.8);
Mass_Ox=OxFR*Mass_Prop/(OxFR+1);  
Mass_Fuel=Mass_Prop-Mass_Ox;
Vol_Ox=Mass_Ox/D_Ox_L;  %Volumen del Oxidante
Vol_Fuel=Mass_Fuel/D_Fuel_L;  %Volumen de combustible
Vol_Ox_Tank=Vol_Ox*1.1;       %Volumen total del tanque, aproximacion
Vol_Fuel_Tank=Vol_Fuel*1.1;   %Volumen total del tanque, aproximacion
Length_Ox_Mass=Vol_Ox/(3.14*.25*Diam_Ox_Tank^2);
Length_Fuel_Mass=Vol_Fuel/(3.14*.25*Diam_Fuel_Tank^2);

Mass_Press_Ox_u=Vol_Ox_Tank*(D_Press_Ox_u)
Mass_Ox_Ox_u=Vol_Ox_Tank*(D_Ox_G_u)

Mass_Press_Fuel_u=Vol_Fuel_Tank*(D_Press_Fuel)
Mass_Fuel_Fuel_u =Vol_Fuel_Tank*(D_Fuel_G)

Vol_Press_Tank=(Mass_Press_Ox_u+Mass_Press_Fuel_u)/(D_Press_Tank_o-D_Press_Tank_u)
Mass_Press_Press_u=Vol_Press_Tank*D_Press_Tank_u
Mass_Press_Press_o=Vol_Press_Tank*D_Press_Tank_o

%% Calcular parametro generales finales

%% Ecuaciones de Flujo

%% Simulacion de la Primera presurizacion al tanque de oxidante
i=1;
dt=0.001;
t(1)=0;
Vol_Ullage_Ox=Vol_Ox_Tank-Vol_Ox;
Mass_Ullage_Ox=Vol_Ullage_Ox*D_Ox_G_o;

Vol_Ullage_Fuel=Vol_Fuel_Tank-Vol_Fuel;
Mass_Ullage_Fuel=Vol_Ullage_Fuel*D_Fuel_G;

Cd_Press_Reg=0.8;
A_Press_Reg=3.14*.25*(0.0025^2);

P_Tank_Ox_t(1)=P_Ox_Vap;
P_Press_Press_t(1)=P_Press_o;
T_Press_Press_t(1)=T_Press_o;
dM_Press_Fuel(1)=0;
dM_Press_Ox(1)=0;
N_Press_Fuel(1)=0;
N_Fuel_Fuel(1)=Mass_Ullage_Fuel/mol_Fuel;

M_Press_Flow(1)=0;
D_Press_Press_t(1)=D_Press_Tank_o;
P_Tank_Fuel_t(1)=P_Fuel_Vap;
P_Fuel_Fuel_t(1)=P_Fuel_Vap;
P_Press_Fuel_t(1)=0;
T_mix_Ox(1)=T_Ox_Tank_o;
T_mix_Fuel(1)=T_Fuel_Tank_o;
T_Press_Fuel_t(1)=T_Press_o*(P_Fuel_Vap/P_Press_o)^((Gamma_Press-1)/Gamma_Press);


D_Press_Ox_mix(1)=0;
D_Ox_mix(1)=D_Ox_G_o;
D_Press_Fuel_mix(1)=0;
D_Fuel_mix(1)=D_Fuel_G;

Mass_Ullage_Ox_t(1)=Mass_Ullage_Ox;
Mass_Ullage_Fuel_t(1)=Mass_Ullage_Fuel;
%P_Tank_Ox_t(i)<P_Tank_Ox
while(i<500)
    %Modelado de un Regulador de Presion
    %Ecuaciones de balance de Masa
    M_Press_Flow(i+1)=Cd_Press_Reg*A_Press_Reg*Pressure_Reg(P_Press_Press_t(i),P_Tank_Fuel_t(i),4*10^6,  T_Press_Press_t(i), Gamma_Press,R_Press);
    dM_Press_Fuel(i+1)=dM_Press_Fuel(i)+M_Press_Flow(i+1)*dt;
    %Mass_Ullage_Ox_t(i+1)=Mass_Ullage_Ox;
    Mass_Ullage_Fuel_t(i+1)=Mass_Ullage_Fuel_t(i);
    
    D_Press_Press_t(i+1)=(Mass_Press_Press_o-dM_Press_Fuel(i+1))/Vol_Press_Tank;
    %D_Press_Fuel_mix(i+1)=(dM_Press_Fuel(i+1))/Vol_Ullage_Fuel;

    N_Press_Fuel(i+1)=dM_Press_Fuel(i+1)/mol_Press;
    N_Fuel_Fuel(i+1)=Mass_Ullage_Fuel_t(i)/mol_Fuel;
    %Ecuaciones de balance de Energia
    u_mol_Mix_Fuel(i)=(hmolar_Press_Press_o*N_Press_Fuel(i+1) + umolar_Fuel_o*N_Fuel_Fuel(i+1))/(N_Press_Fuel(i+1)+N_Fuel_Fuel(i+1));
    %[T_mix_Fuel(i+1),D_Fuel_mix(i+1),D_Press_Fuel_mix(i+1)]=Mix_Equilibrium_1(H_Mix_Fuel(i),Mass_Ullage_Fuel_t(i),dM_Press_Fuel(i+1),Vol_Ullage_Fuel,Fuel,Pressurant,D_Fuel_mix(i),T_mix_Fuel(i));  
    
    P_Press_Fuel_t(i+1)=(((hmolar_Press_Press_o*N_Press_Fuel(i+1) + (umolar_Fuel_o-u_mol_Mix_Fuel(i))*N_Fuel_Fuel(i+1))*(R/(Vol_Ullage_Fuel*u_mol_Mix_Fuel(i)))*T_Press_Press_t(i)*(P_Press_Press_t(i))^(-1+1/Gamma_Press)))^Gamma_Press;
    P_Tank_Fuel_t(i+1)=P_Press_Fuel_t(i+1)+P_Fuel_Fuel_t(i);
    T_Press_Fuel_t(i+1)=T_Press_o*(P_Tank_Fuel_t(i+1)/P_Press_o)^((Gamma_Press-1)/Gamma_Press);
    T_mix_Fuel(i+1)=P_Tank_Fuel_t(i+1)*Vol_Ullage_Fuel/((N_Press_Fuel(i+1)+N_Fuel_Fuel(i+1))*R);
    %P_Press_Fuel_t(i+1)=P_Press_Fuel_t(i+1)*T_mix_Fuel(i+1)/T_Press_Fuel_t(i+1);
    P_Fuel_Fuel_t(i+1) =P_Fuel_Vap*T_mix_Fuel(i+1)/T_Fuel_Tank_o;
    
    % Nuevos Estados
    %P_Press_Fuel_t(i+1)=py.CoolProp.CoolProp.PropsSI('P','H',h_Press_Press_o,  'D',D_Press_Fuel_mix(i+1),  Pressurant);
    %T_Press_Fuel_t(i+1)=py.CoolProp.CoolProp.PropsSI('T','H',h_Press_Press_o,  'D',D_Press_Fuel_mix(i+1),  Pressurant);

    %py.CoolProp.CoolProp.PropsSI('P','T',T_mix_Fuel(i+1),  'D',D_Fuel_mix(i+1),  Fuel);

    P_Press_Press_t(i+1)=py.CoolProp.CoolProp.PropsSI('P','H',h_Press_Press_o,  'D',D_Press_Press_t(i+1),  Pressurant);
    T_Press_Press_t(i+1)=py.CoolProp.CoolProp.PropsSI('T','H',h_Press_Press_o,  'D',D_Press_Press_t(i+1),  Pressurant);
    t(i+1)=t(i)+dt;

    i=i+1;
end
close all
figure(1)
plot(t,P_Press_Press_t, t, P_Fuel_Fuel_t, t, P_Press_Fuel_t, t, P_Tank_Fuel_t)
legend('P Press Press', 'P Fuel Fuel', 'P Press Fuel', 'P Fuel Tank')
title('Presiones')

figure(2)
plot(t,T_Press_Press_t ,t, T_Press_Fuel_t, t, T_mix_Fuel )
legend('T Press Press', 'T Press Fuel', 'T Mix Fuel')
title('Temperaturas')

figure(3)
plot(t,dM_Press_Fuel)
legend('Press Mass Gain')
title('Masas')

figure(4)
plot(t,M_Press_Flow)
legend('Press Mass Flow')
title('Flujo de Masas')



