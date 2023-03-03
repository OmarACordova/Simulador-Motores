%% Hacer Simulacion
clear all
close all

syms b1_Tvap b2_Tvap b3_Tvap b4_Tvap Pcrita Tr2
b1_Tvap=-6.71893;
b2_Tvap=+1.35966;
b3_Tvap=-1.3779;
b4_Tvap=-4.051;
Pcrita=7.25*10^6;

b1_hlsat=-200;
b2_hlsat=+116.043
b3_hlsat=-917.225
b4_hlsat=+794.779
b5_hlsat=-589.587

b1_hgsat=-200
b2_hgsat=+440.055
b3_hgsat=-459.701
b4_hgsat=+434.081
b5_hgsat=-485.338

b1_Cpg=+132.632
b2_Cpg=+0.052187
b3_Cpg=-0.364923
b4_Cpg=-1.20233
b5_Cpg=+0.536141

b1_Cpl=+2.49973
b2_Cpl=+0.023454
b3_Cpl=-3.80136
b4_Cpl=+13.0945
b5_Cpl=-14.5180


Alg=((6*0.0254)^2)*0.25*3.14;
Ta1=10+273
Tl1=-15+273
Tra=Ta1/(273+36)
Trl=Tl1/(273+36)
Pa1=3.127*10^6;
hlg=100;
Vg1=20*10^-3;
Vl1=1*10^-3;

mma=0.044;
R=8.314;
Rola=907;
dt=0.75;
time=0;

eqPTsaturation=Pa1==Pcrita*exp((b1_Tvap*(1-Tr2)+b2_Tvap*(1-Tr2).^(1.5)+b3_Tvap*(1-Tr2).^(2.5)+b4_Tvap*(1-Tr2).^(5))./Tr2);
Tvapa2=vpasolve(eqPTsaturation, Tr2)*(273+36)
Trvap=Tvapa2/(273+36);

ugsat1=b1_hgsat + b2_hgsat*(1-Trvap)^(1/3) + b3_hgsat*(1-Trvap)^(2/3) + b4_hgsat*(1-Trvap) + b5_hgsat*(1-Trvap)^(4/3)
ulsat1=b1_hlsat + b2_hlsat*(1-Trvap)^(1/3) + b3_hlsat*(1-Trvap)^(2/3) + b4_hlsat*(1-Trvap) + b5_hlsat*(1-Trvap)^(4/3)
hvapa=ugsat1-ulsat1;
Cvla=b1_Cpl*(1 + b2_Cpl*(1-Trvap)^-1 + b3_Cpl*(1-Trvap) + b4_Cpl*(1-Trvap)^2 + b5_Cpl*(1-Trvap)^3)
Cvga=b1_Cpg*(1 + b2_Cpg*(1-Trvap)^(-2/3) + b3_Cpg*(1-Trvap)^(-1/3) + b4_Cpg*(1-Trvap)^(1/3) + b5_Cpg*(1-Trvap)^(2/3))
uga1=ugsat1+Cvla*(Ta1-Tvapa2)
ula1=ulsat1+Cvga*(Tl1-Tvapa2)

nga1 =(Pa1*Vg1)/(R*Ta1)
nla1 =(Rola*Vl1)/mma

Pvapa1=Pcrita*exp((b1_Tvap*(1-Trl)+b2_Tvap*(1-Trl).^(1.5)+b3_Tvap*(1-Trl).^(2.5)+b4_Tvap*(1-Trl).^(5))./Trl);


Vtank =Vg1 + Vl1

dqlg=dt*Alg*hlg*(Ta1-Tl1)

datos(1,:)=[time, nga1, nla1, uga1, ula1, dqlg, Ta1, Tvapa2, Tl1, Pa1, Pvapa1, Vg1, Vl1]
i=1;

%%

if(Ta1>Tvapa2 && Tl1<Tvapa2 )
    while(time<75)
    
        nga2 =(hvapa*nga1 - Alg*Ta1*hlg + Alg*Tl1*hlg + Cvga*Ta1*nga1 - Cvla*Tl1*nga1 - Cvga*Tvapa2*nga1 + Cvla*Tvapa2*nga1)/(hvapa + Cvga*Ta1 - Cvla*Tl1 - Cvga*Tvapa2 + Cvla*Tvapa2)
        nla2 =(hvapa*nla1 + Alg*Ta1*hlg - Alg*Tl1*hlg + Cvga*Ta1*nla1 - Cvla*Tl1*nla1 - Cvga*Tvapa2*nla1 + Cvla*Tvapa2*nla1)/(hvapa + Cvga*Ta1 - Cvla*Tl1 - Cvga*Tvapa2 + Cvla*Tvapa2)
        Ta2  =(Cvga*Cvla*Ta1^2*nga1 - Alg*Cvla*Ta1^2*hlg - Cvla^2*Ta1*Tl1*nga1 + Cvla^2*Ta1*Tvapa2*nga1 + Cvla*Ta1*hvapa*nga1 + Alg*Ta1*hlg*uga1 - Alg*Tl1*hlg*uga1 - Alg*Ta1*dt*hlg*hvapa + Alg*Tl1*dt*hlg*hvapa - Alg*Cvga*Ta1^2*dt*hlg - Alg*Cvla*Tl1^2*dt*hlg + Alg*Cvla*Ta1*Tl1*hlg - Cvga*Cvla*Ta1*Tvapa2*nga1 + Alg*Cvga*Ta1*Tl1*dt*hlg + Alg*Cvla*Ta1*Tl1*dt*hlg + Alg*Cvga*Ta1*Tvapa2*dt*hlg - Alg*Cvla*Ta1*Tvapa2*dt*hlg - Alg*Cvga*Tl1*Tvapa2*dt*hlg + Alg*Cvla*Tl1*Tvapa2*dt*hlg)/(Cvla*(hvapa*nga1 - Alg*Ta1*hlg + Alg*Tl1*hlg + Cvga*Ta1*nga1 - Cvla*Tl1*nga1 - Cvga*Tvapa2*nga1 + Cvla*Tvapa2*nga1))
        Tl2  =(Cvga^2*Ta1*Tl1*nla1 - Cvga*Cvla*Tl1^2*nla1 - Alg*Cvga*Tl1^2*hlg - Cvga^2*Tl1*Tvapa2*nla1 + Cvga*Tl1*hvapa*nla1 - Alg*Ta1*hlg*ula1 + Alg*Tl1*hlg*ula1 + Alg*Ta1*dt*hlg*hvapa - Alg*Tl1*dt*hlg*hvapa + Alg*Cvga*Ta1^2*dt*hlg + Alg*Cvla*Tl1^2*dt*hlg + Alg*Cvga*Ta1*Tl1*hlg + Cvga*Cvla*Tl1*Tvapa2*nla1 - Alg*Cvga*Ta1*Tl1*dt*hlg - Alg*Cvla*Ta1*Tl1*dt*hlg - Alg*Cvga*Ta1*Tvapa2*dt*hlg + Alg*Cvla*Ta1*Tvapa2*dt*hlg + Alg*Cvga*Tl1*Tvapa2*dt*hlg - Alg*Cvla*Tl1*Tvapa2*dt*hlg)/(Cvga*(hvapa*nla1 + Alg*Ta1*hlg - Alg*Tl1*hlg + Cvga*Ta1*nla1 - Cvla*Tl1*nla1 - Cvga*Tvapa2*nla1 + Cvla*Tvapa2*nla1))
        dqlg =Alg*dt*hlg*(Ta2 - Tl2)
    
        Vl2=nla2*mma/Rola;
        Vg2=Vg1+Vl1-Vl2;
    
        %Actualizar Variables
        nga1=nga2;
        nla1=nla2;
        ula1=ula1+Cvla*(Tl2-Tl1);
        uga1=uga1+Cvla*(Ta2-Ta1);
        Ta1=Ta2;
        Pa2=R*Ta2*nga2/Vg2;
        Vg1=Vg2;
        Vl1=Vl2;

        eqPTsaturation=Pa2==Pcrita*exp((b1_Tvap*(1-Tr2)+b2_Tvap*(1-Tr2).^(1.5)+b3_Tvap*(1-Tr2).^(2.5)+b4_Tvap*(1-Tr2).^(5))./Tr2);
        Tvapa2=vpasolve(eqPTsaturation, Tr2)*(273+36);
        Trvap=Tvapa2/(273+36);
        
        ugsat1=b1_hgsat + b2_hgsat*(1-Trvap)^(1/3) + b3_hgsat*(1-Trvap)^(2/3) + b4_hgsat*(1-Trvap) + b5_hgsat*(1-Trvap)^(4/3)
        ulsat1=b1_hlsat + b2_hlsat*(1-Trvap)^(1/3) + b3_hlsat*(1-Trvap)^(2/3) + b4_hlsat*(1-Trvap) + b5_hlsat*(1-Trvap)^(4/3)
        hvapa=ugsat1-ulsat1;
        Cvla=b1_Cpl*(1 + b2_Cpl*(1-Trvap)^-1 + b3_Cpl*(1-Trvap) + b4_Cpl*(1-Trvap)^2 + b5_Cpl*(1-Trvap)^3)
        Cvga=b1_Cpg*(1 + b2_Cpg*(1-Trvap)^(-2/3) + b3_Cpg*(1-Trvap)^(-1/3) + b4_Cpg*(1-Trvap)^(1/3) + b5_Cpg*(1-Trvap)^(2/3))
        Trl=Tl2/(273+36);
        Pvapa1=Pa2==Pcrita*exp((b1_Tvap*(1-Trl)+b2_Tvap*(1-Trl).^(1.5)+b3_Tvap*(1-Trl).^(2.5)+b4_Tvap*(1-Trl).^(5))./Trl);
        uga1=ugsat1+Cvla*(Ta1-Tvapa2)
        ula1=ulsat1+Cvga*(Tl1-Tvapa2)
    
        
    
        time=time+dt
        i=i+1;
    
        datos(i,:)=[time, nga2, nla2, uga1, ula1, dqlg, Ta2, Tvapa2, Tl2,Pa2, Pvapa1, Vg2, Vl2];
    end
%%
    figure(1)
    plot(datos(:,1), datos(:,7), datos(:,1), datos(:,8),datos(:,1), datos(:,9))
    legend('Temp Vapor', 'Temp Boil', 'Temp Liq')
    
    figure(2)
    plot(datos(:,1), datos(:,2), datos(:,1), datos(:,3))
    legend('Mol Vapor', 'Mol Liq')
    
    figure(3)
    plot(datos(:,1), datos(:,6))
    legend('Heat Flow')
    
    figure(4)
    plot(datos(:,1), datos(:,10), datos(:,1), datos(:,11))
    legend('Pres A', 'Pres Vap Liq')
    
    figure(5)
    plot(datos(:,1), datos(:,12), datos(:,1), datos(:,13))
    legend('Vol Vapor', 'Vol Liq')

else
    Mensaje='Error: COmbiancion de Ta, Tl y Pa son imposibles'

end

