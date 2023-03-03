close all
clear all
b1= 132.632;
b2= 0.052187;
b3=-0.364923
b4= -1.20233
b5= 0.53614
T=[-90:20]
Tr=(273+T)/(36+273)
Cpgsat=b1*[1+b2*(1-Tr).^(-2/3)+b3*(1-Tr).^(-1/3)+b4*(1-Tr).^(1/3)+b5*(1-Tr).^(2/3)]

hvapa=230;
mma=0.044;
R=8.314;



%% 
clear all
close all
syms nga1 nga2 nla1 nla2              %Moles de gases y liquidos 
syms uga1 ula1 uga2 ula2 Cvla Cvga    %Energias Internas y Entalpias 
syms Qlg1 Qlg2 dqlg hvapa Alg                    %Tranferencia de Calor Vapor liquido
syms Ta1 Ta2 Tvapa1 Tvapa2 Tl1 Tl2    %Temperaturas de Transferencia de calor Vapor-Liquido
syms hlg                              %Coeficientes de conveccion
syms Pa1 Pa2 R Vtank Vg1 Vg2 Vl1 Vl2 mma
syms Rola
%syms Pb0 Tb0 Gamma
syms dt
syms Tvaparef Pref

%nga2, Ta2, dqlg, nla2, Tl2, Vl2, 
eq1=nga2*(uga1+Cvga*(Ta2-Ta1))==(nga1*uga1-dqlg)           
eq2=nla2*(ula1+Cvla*(Tl2-Tl1))==(nla1*ula1+dqlg)

eq3=dqlg==dt*(nla2-nla1)*(Cvga*(Ta2-Tvapa2)+hvapa+Cvla*(Tvapa2-Tl2))
eq4=dqlg==dt*Alg*hlg*(Ta2-Tl2)
%conservacion de al masa
eq5=nga1+nla1==nga2+nla2

%Convervacion del Volumen
eq6=Vl1+Vg1==Vl2+Vg2

EQ=[eq1,eq2,eq3,eq4,eq5,eq6]
X=[nga2, Ta2, dqlg, nla2, Tl2, Vl2]
Sol=solve(EQ,X)


eqQ=Alg*hlg*(Ta2-Tl2)==(nla2-nla1)*(Cvga*(Ta2-Tvapa2)+hvapa+Cvla*(Tvapa2-Tl2))

Tl2=((nla1 - nla2)*(hvapa + Cvga*(Ta2 - Tvapa2) + Cvla*Tvapa2) + Alg*Ta2*hlg)/(Alg*hlg + Cvla*(nla1 - nla2))
%%
syms b1 b2 b3 b4 Pvapa1 Pcrita Tr1
b1=-6.71893;
b2=+1.35966;
b3=-1.3779;
b4=-4.051;
Pcrita=7.25*10^6;
Pvapa1=4*10^6;

nga2 =(hvapa*nga1 - Alg*Ta1*hlg + Alg*Tl1*hlg + Cvga*Ta1*nga1 - Cvla*Tl1*nga1 - Cvga*Tvapa2*nga1 + Cvla*Tvapa2*nga1)/(hvapa + Cvga*Ta1 - Cvla*Tl1 - Cvga*Tvapa2 + Cvla*Tvapa2)
nla2 =(hvapa*nla1 + Alg*Ta1*hlg - Alg*Tl1*hlg + Cvga*Ta1*nla1 - Cvla*Tl1*nla1 - Cvga*Tvapa2*nla1 + Cvla*Tvapa2*nla1)/(hvapa + Cvga*Ta1 - Cvla*Tl1 - Cvga*Tvapa2 + Cvla*Tvapa2)
Ta2 =(Cvga^2*Ta1^2*nla1 - hvapa*nla1*uga1 + hvapa*nla1*ula1 + Alg*Cvga*Ta1^2*hlg - Cvga^2*Ta1*Tvapa2*nla1 + Cvga*Ta1*hvapa*nla1 - Alg*Ta1*hlg*uga1 + Alg*Tl1*hlg*uga1 - Cvga*Ta1*nla1*uga1 + Cvga*Ta1*nla1*ula1 + Cvla*Tl1*nla1*uga1 - Cvla*Tl1*nla1*ula1 + Cvga*Tvapa2*nla1*uga1 - Cvla*Tvapa2*nla1*uga1 - Cvga*Tvapa2*nla1*ula1 + Cvla*Tvapa2*nla1*ula1 + Alg*Ta1*dt*hlg*hvapa - Alg*Tl1*dt*hlg*hvapa + Alg*Cvga*Ta1^2*dt*hlg + Alg*Cvla*Tl1^2*dt*hlg - Alg*Cvga*Ta1*Tl1*hlg - Cvga*Cvla*Ta1*Tl1*nla1 + Cvga*Cvla*Ta1*Tvapa2*nla1 - Alg*Cvga*Ta1*Tl1*dt*hlg - Alg*Cvla*Ta1*Tl1*dt*hlg - Alg*Cvga*Ta1*Tvapa2*dt*hlg + Alg*Cvla*Ta1*Tvapa2*dt*hlg + Alg*Cvga*Tl1*Tvapa2*dt*hlg - Alg*Cvla*Tl1*Tvapa2*dt*hlg)/(Cvga*(hvapa*nla1 + Alg*Ta1*hlg - Alg*Tl1*hlg + Cvga*Ta1*nla1 - Cvla*Tl1*nla1 - Cvga*Tvapa2*nla1 + Cvla*Tvapa2*nla1))
Tl2 =(hvapa*nga1*uga1 - Cvla^2*Tl1^2*nga1 - hvapa*nga1*ula1 + Alg*Cvla*Tl1^2*hlg + Cvla^2*Tl1*Tvapa2*nga1 + Cvla*Tl1*hvapa*nga1 + Alg*Ta1*hlg*ula1 - Alg*Tl1*hlg*ula1 + Cvga*Ta1*nga1*uga1 - Cvga*Ta1*nga1*ula1 - Cvla*Tl1*nga1*uga1 + Cvla*Tl1*nga1*ula1 - Cvga*Tvapa2*nga1*uga1 + Cvla*Tvapa2*nga1*uga1 + Cvga*Tvapa2*nga1*ula1 - Cvla*Tvapa2*nga1*ula1 - Alg*Ta1*dt*hlg*hvapa + Alg*Tl1*dt*hlg*hvapa - Alg*Cvga*Ta1^2*dt*hlg - Alg*Cvla*Tl1^2*dt*hlg - Alg*Cvla*Ta1*Tl1*hlg + Cvga*Cvla*Ta1*Tl1*nga1 - Cvga*Cvla*Tl1*Tvapa2*nga1 + Alg*Cvga*Ta1*Tl1*dt*hlg + Alg*Cvla*Ta1*Tl1*dt*hlg + Alg*Cvga*Ta1*Tvapa2*dt*hlg - Alg*Cvla*Ta1*Tvapa2*dt*hlg - Alg*Cvga*Tl1*Tvapa2*dt*hlg + Alg*Cvla*Tl1*Tvapa2*dt*hlg)/(Cvla*(hvapa*nga1 - Alg*Ta1*hlg + Alg*Tl1*hlg + Cvga*Ta1*nga1 - Cvla*Tl1*nga1 - Cvga*Tvapa2*nga1 + Cvla*Tvapa2*nga1))
dqlg =Alg*dt*hlg*(Ta2 - Tl2)
Vl2 =Vg1 - Vg2 + Vl1

nga1=nga2
nla1=nla2
ula1=ula1+Cvla*(Tl2-Tl1)
uga1=uga1+Cvla*(Ta2-Ta1)
Ta1=Ta2
Tr2=Ta2/(273+36)

eqPTsaturation=Pvapa1==Pcrita*exp((b1*(1-Tr2)+b2*(1-Tr2).^(1.5)+b3*(1-Tr2).^(2.5)+b4*(1-Tr2).^(5))./Tr2);
Tvapa2=vpasolve(eqPTsaturation, Tr1)*(273+36);

%%
b1_Tvap=-6.71893;
b2_Tvap=+1.35966;
b3_Tvap=-1.3779;
b4_Tvap=-4.051;
Pcrita=7.25*10^6;

a_VW_NO2=0.5354;  % Atraccion entre moelesuclas de Van der Walls
b_VW_NO2=0.00004424;  % Tamaño de las Moleculas Van der Walls
a_VW_N2=0.139
b_VW_N2=0.0000391
n=linspace(-10,30,100)
Vg=5*10^-3
Tg=200
R=8.314

Trl=Tg/(273+36)

Pvapl=Pcrita*exp((b1_Tvap*(1-Trl)+b2_Tvap*(1-Trl).^(1.5)+b3_Tvap*(1-Trl).^(2.5)+b4_Tvap*(1-Trl).^(5))./Trl);

P=n*R*Tg./(Vg-n*b_VW_N2)-a_VW_N2*(n/Vg).^2
P2=n*R*Tg/(Vg)
plot(n, P, n, P2,[-10,20],[Pvapl, Pvapl] )



















