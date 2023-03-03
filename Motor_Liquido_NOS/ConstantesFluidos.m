
% Temperatura del Liquido en funcion de la Presion parcial de vapor
% Tvap=b1_TPvap*exp(b2_TPvap*PNos1) +  b3_TPvap*exp(-5.139e-07*PNos1)
% global b1_TPvap b2_TPvap b3_TPvap b4_TPvap

% ---- NOS -----
b1_TPvap=268.2;
b2_TPvap=2.067e-08;
b3_TPvap=-64.09;
b4_TPvap=-5.139e-07;

% ---- CO2 ------
% b1_TPvap=258.8;
% b2_TPvap=2.237e-08;
% b3_TPvap= -61.55;
% b4_TPvap=-6.496e-07;


% Presion de vapor en funcion de la Temperatura del liquido
% PNos1=b1_PvapT*Te1^3 + b2_PvapT*Te1^2 + b3_PvapT*Te1 + b4_PvapT;
% global b1_PvapT b2_PvapT b3_PvapT b4_PvapT

% ---- NOS-------
b1_PvapT=2.792;
b2_PvapT=-1527; 
b3_PvapT=2.898e+05;
b4_PvapT=-1.899e+07;

% -----CO2------
% b1_PvapT=3.922;
% b2_PvapT=-2283; 
% b3_PvapT=4.619e+05;
% b4_PvapT=-3.231e+07;

% Entalpia/Energia Interna del oxidante Liquido en funcion de la temperatura
% global b1_hlsat b2_hlsat b3_hlsat b4_hlsat b5_hlsat
% ulsat1=1000*(b1_hlsat + b2_hlsat*(1-Trl)^(1/3) + b3_hlsat*(1-Trl)^(2/3) + b4_hlsat*(1-Trl) + b5_hlsat*(1-Trl)^(4/3));   % cambiar

% ---- NOS ----
b1_hlsat=-200;
b2_hlsat=+116.043;
b3_hlsat=-917.225;
b4_hlsat=+794.779;
b5_hlsat=-589.587;

% ---- CO2 ----
% b1_hlsat= 318.4;
% b2_hlsat=-52.68;
% b3_hlsat=-416.4;
% b4_hlsat= 77.48;
% b5_hlsat=-233.3;

% Entalpia/Energia Interna del oxidante gaseoso en funcion de la temperatura
% ugsat1=1000*(b1_hgsat + b2_hgsat*(1-Trg)^(1/3) + b3_hgsat*(1-Trg)^(2/3) + b4_hgsat*(1-Trg) + b5_hgsat*(1-Trg)^(4/3));   % cambiar
% global b1_hgsat b2_hgsat b3_hgsat b4_hgsat b5_hgsat
% ---- NOS ----
b1_hgsat=-200;
b2_hgsat=+440.055;
b3_hgsat=-459.701;
b4_hgsat=+434.081;
b5_hgsat=-485.338;

% ---- CO2 ----
% b1_hgsat=347.9;
% b2_hgsat=82.86;
% b3_hgsat=520.1;
% b4_hgsat=-646;
% b5_hgsat=-69.64;

% Entalpia del oxidante liquido en funcion de la presion parcial de vapor
% Hl0 = (b1_HlP)*PNos1.^3 + (b2_HlP)*PNos1.^2 + (b3_HlP).*PNos1 + b4_HlP;                     % Sustancia
% global b1_HlP b2_HlP b3_HlP b4_HlP

% --- NOS -----
b1_HlP=5.842e-16;
b2_HlP=-8.466e-09;
b3_HlP=0.06278;
b4_HlP=-4.315e+05;

% --- CO2 ----
% b1_HlP=9.387e-16;
% b2_HlP=-1.205e-8;
% b3_HlP=7.435e-02 ;
% b4_HlP=48.02e+3;

% Entalpia del oxidante gas en funcion de la presion
% parcial de vapor
% Hg0 = (b1_HgP)*Pc^8 + (b2_HgP)*Pc^7 + (b3_HgP)*Pc^6 + (b4_HgP)*Pc^5 + ...
% (b5_HgP)*Pc^4 + (b6_HgP)*Pc^3 + (b7_HgP)*Pc^2 + (b8_HgP)*Pc + (b9_HgP);  
% global b1_HgP b2_HgP b3_HgP b4_HgP b5_HgP b6_HgP b7_HgP b8_HgP b9_HgP

% ---- NOS -----
b1_HgP=-1.223e-47; 
b2_HgP=3.658e-40; 
b3_HgP=-4.626e-33; 
b4_HgP=3.221e-26; 
b5_HgP=-1.346e-19;
b6_HgP=3.452e-13;
b7_HgP=-5.306e-07; 
b8_HgP=0.4495;
b9_HgP=-2.311e+05;

% ---- CO2 ------
% b1_HgP=-1.532e-48; 
% b2_HgP= 4.282e-41; 
% b3_HgP=-5.041e-34; 
% b4_HgP= 3.26e-27; 
% b5_HgP=-1.272e-20;
% b6_HgP= 3.107e-14;
% b7_HgP=-4.933e-8; 
% b8_HgP= 4.845e-02;
% b9_HgP= 415e+3 ;

% Capacidad termica del oxidante gas en Funcion de la temperatura
% Cpg0 = b1_Cpg*exp(b2_Cpg*TOx) +  (b3_Cpg)*exp(b4_Cpg *TOx);
% global b1_Cpg b2_Cpg b3_Cpg b4_Cpg

% --- NOS -----
b1_Cpg=308.4;
b2_Cpg=0.005439; 
b3_Cpg=6.942e-13; 
b4_Cpg=0.121;

% ---- CO2 ----
% b1_Cpg=333;
% b2_Cpg=0.003165; 
% b3_Cpg=5.948e-13 ; 
% b4_Cpg=  0.1155;

% Capacidad termica del oxidante liquido en Funcion de la temperatura
% Cpl0 = b1_Cpl*exp( b2_Cpl *TOx) + (b3_Cpl)*exp(b4_Cpl  *TOx);                                            % Sustancia
% global b1_Cpl b2_Cpl b3_Cpl b4_Cpl
% --- NOS ----
b1_Cpl=665;
b2_Cpl=0.004254;
b3_Cpl=8.183e-10;
b4_Cpl=0.09463;

% ---- CO2 -----
% b1_Cpl= 4880;
% b2_Cpl=-0.003656;
% b3_Cpl=3.719e-11;
% b4_Cpl=0.1097;

% Densidad del oxidante liquido en funcion de la temperatura
% RolOx=RolOx_c*exp( b1_RolOx*(1-Trl)^(1/3) + b2_RolOx*(1-Trl)^(2/3) + b3_RolOx*(1-Trl) + b4_RolOx*(1-Trl)^(4/3) );
% global b1_RolOx b2_RolOx b3_RolOx b4_RolOx

% ---- NOS ----
b1_RolOx=+1.72328;
b2_RolOx=-0.83950;
b3_RolOx=+0.51060;
b4_RolOx=-0.10412;

% ---- CO2 ----
% b1_RolOx= 1.753;
% b2_RolOx=-0.8536;
% b3_RolOx= 0.6431;
% b4_RolOx=-0.2339;

% Volumen Especifico de Oxidante Liquido en funcion de la presion parcial
% de vapor del oxidante
% Vl0 = b1_VlP*exp((b2_VlP)*PNos1) +  (b3_VlP)*exp((b4_VlP)*PNos1);
% global b1_VlP b2_VlP b3_VlP b4_VlP 
% --- NOS ---
b1_VlP=0.0008209;
b2_VlP=7.749e-08;
b3_VlP=2.269e-09;
b4_VlP=1.649e-06;

% --- CO2 ---
%b1_VlP=0.0008359;
%b2_VlP=7.327e-08;
%b3_VlP=1.852e-10;
%b4_VlP=1.999e-06;

% Volumen Especifico del Oxdante Gas en funcion de la presion parcial
% de vapor del oxidante
% global b1_VgP b2_VgP b3_VgP b4_VgP
% Vg2I =  b1_VgP*exp((b2_VgP )*PcI) +  (b3_VgP)*exp((b4_VgP)*PcI);

% ---- NOS ----
b1_VgP= 0.1382 ;
b2_VgP=-2.325e-06; 
b3_VgP=0.03549;
b4_VgP=-3.497e-07 ;

% ---- CO2 ----
%b1_VgP=0.03712;
%b2_VgP=-1.501e-06 ; 
%b3_VgP=0.03537;
%b4_VgP=-3.576e-07;

a_VW_N2=0.139;          %Paramtros Van der Waals Nitrogeno
b_VW_N2=0.0000391;     %Paramtros van der Waals Nitrogeno 

Tatm=273+25;
Patm=77*10^3;
g=9.8;
R=8.314;

Gammag=1.4;
mmg=0.028;
Kg=26*10^-3;
Cvg=740;
Mug=18*10^-6;
Rg=R/mmg;

mmv=0.044;
RolOx_c=467;
Beta_lOx=1/128;
Muv=1.6832e-05	;
Mul_Ox= 7.3870e-05	;
KlOx=0.090130;
Kv=0.029961;
TfMin=243;
Pcrita=7.25*10^6;
Tcrita=273+31;
Rv=R/mmv;

RolFu=790;              % Densidad del Etanol
Beta_lFu=1/912;
CvFu=118;
Mul_Fu=  0.00110;
KlFu=167*10^-3;
% .---------------------------------------------------------------------
L_char=1.69;            %Sustancia

% (To_1*OFR(i).^2 + To_2*OFR(i) + To_3)./(OFR(i).^3 + To_4*OFR(i).^2 + To_5*OFR(i) + To_6)
To_1=-4.765e+10;
To_2= 1.042e+11;
To_3=-1.742e+11;
To_4=-1.808e+07;
To_5= 6.731e+07;
To_6=-1.464e+08;

% (Mm_1*OFR(i).^2 +    Mm_2.*OFR(i) + Mm_3)./ (OFR(i).^2 + Mm_4.*OFR(i) + Mm_5);
Mm_1=0.0292;
Mm_2=-0.06132;
Mm_3=0.09785;
Mm_4=-1.92;
Mm_5=4.573;

% (Gamma_1.*OFR(i).^3 + Gamma_2.*OFR(i).^2 + Gamma_3.*OFR(i) + Gamma_4)./(OFR(i).^2 + Gamma_5.*OFR(i) + Gamma_6);
Gamma_1=0.004028;
Gamma_2=1.185;
Gamma_3=-2.961;
Gamma_4=2.672;
Gamma_5=-2.527;
Gamma_6=2.241;

% ----------------------------------------------------------------------
Cvw=500;
Ktank=15;
Xtank=10*10^-3;
% ---------------------------------------------------------------------
Al6063_Ro=2700;
Al6063_Sigma=250*10^6;

Al6061_Ro=2700;
Al6061_Sigma=280*10^6;

SS304_Ro=7900;
SS304_Sigma=250*10^6;

Cu_Sigma=90*10^6;
Cu_Ro=9400;

CFRP_Sigma=1000*10^6;
CFRP_Ro=2000;
