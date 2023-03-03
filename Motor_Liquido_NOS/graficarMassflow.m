clear all
close all
TOx=273+5
TlOx=TOx;
RolOx_c=452;
Cd=0.66
b1_RolOx=+1.72328;
b2_RolOx=-0.83950;
b3_RolOx=+0.51060;
b4_RolOx=-0.10412;

b1_Cpl=+2.49973;
b2_Cpl=+0.023454;
b3_Cpl=-3.80136;
b4_Cpl=+13.0945;
b5_Cpl=-14.5180;

LInjOx=10*10^-3;
DInjOx=2*10^-3
NInjOx=15;
AInjOx=NInjOx*3.14*.25*DInjOx^2;

LInjFu=10*10^-3;
DInjFu=2*10^-3
NInjFu=3;
AInjFu=NInjFu*3.14*.25*DInjFu^2

Cd=0.66

WLineOx=0.889*10^-3;
DLineOx=0.0254*1/2-2*WLineOx;     % Variable de Diseño
WLineFu=0.762*10^-3;

DLineFu=0.0254*3/8-2*WLineFu;     % Variable de Diseño
LLineOx=1;
LLineFu=1;
Mul_Ox=  0.000077;
MulFu=1.074^-3;

ALineFu=DLineFu^2*0.25*3.14
kLineFu=3*10^-6/DLineFu;            % Rugosidad relativa de la tuberia de cobre o acero inoxidable flexible "Cooper tubing"
ALineOxo=DLineOx^2*0.25*3.14;
kLineOx=3*10^-6/DLineOx;            % Rugosidad relativa de la tuberia de cobre o acero inoxidable flexible "Cooper tubing"
Trl=TlOx/(273+36);               % Temepratura reducuida de los liquidos del oxidnante
RolFu=790;                          % Desnidad del Etanol

% Densidad liquida del Oxidante en funcion de la temperatura
RolOx=RolOx_c*exp(b1_RolOx*(1-Trl)^(1/3) + b2_RolOx*(1-Trl)^(2/3) + b3_RolOx*(1-Trl) + b4_RolOx*(1-Trl)^(4/3)) 

Cvl=1000*b1_Cpl*(1 + b2_Cpl*(1-Trl)^-1 + b3_Cpl*(1-Trl) + b4_Cpl*(1-Trl)^2 + b5_Cpl*(1-Trl)^3);


Mul_Fu=  0.00110;
Pumbling(1)=DLineOx;
Pumbling(2)=DLineFu;
Pumbling(3)=LLineOx;
Pumbling(4)=LLineFu;
Pumbling(5)=Mul_Ox;
Pumbling(6)=Mul_Fu;
Pumbling(7)=kLineFu;
Pumbling(8)=kLineOx;
Pumbling(9)=Cvl;
Pumbling(10)=RolFu;
Pumbling(11)=Cd;
Pumbling(12)=AInjOx;
Pumbling(13)=AInjFu;

i=1;
j=1;

for PoOx=  1*10^6:1*10^5:6*10^6
    j=1;
    for Pc=1*10^6:5*10^4:6*10^6
        [G_Dyer(i,j) G_Omega(i,j) G_SPI(i,j) G_Omar(i,j)]=MassFlow_InjLine(PoOx, PoOx,Pc, TOx, Pumbling);
        j=j+1;
    end
    i=i+1;
end

[PoOx, Pc]=meshgrid(1*10^6:1*10^5:6*10^6, 1*10^6:5*10^4:6*10^6);

figure(1)
plot3(PoOx,Pc, G_Dyer)

figure(2)
plot3(PoOx,Pc, G_Omega)

figure(3)
plot3(PoOx,Pc, G_SPI)

figure(4)
plot3(PoOx,Pc, G_Omar)


