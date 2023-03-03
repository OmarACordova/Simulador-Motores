

% T=b1_TPvap*exp(b2_TPvap*PNos1) +  b3_TPvap*exp(-5.139e-07*PNos1)
% global b1_TPvap b2_TPvap b3_TPvap b4_TPvap
b1_TPvap=268.2;
b2_TPvap=2.067e-08;
b3_TPvap=-64.09;
b4_TPvap=-5.139e-07;

% PNos1=b1_PvapT*Te1^3 + b2_PvapT*Te1^2 + b3_PvapT*Te1 + b4_PvapT;
% global b1_PvapT b2_PvapT b3_PvapT b4_PvapT
b1_PvapT=2.792;
b2_PvapT=-1527; 
b3_PvapT=2.898e+05;
b4_PvapT=-1.899e+07;

TfMin=223;
Pcrita=7.25*10^6;
Tcrita=273+36;

% global b1_hlsat b2_hlsat b3_hlsat b4_hlsat b5_hlsat
b1_hlsat=-200;
b2_hlsat=+116.043;
b3_hlsat=-917.225;
b4_hlsat=+794.779;
b5_hlsat=-589.587;

% 1000*(b1_hgsat + b2_hgsat*(1-Trg)^(1/3) + b3_hgsat*(1-Trg)^(2/3) + b4_hgsat*(1-Trg) + b5_hgsat*(1-Trg)^(4/3));
% global b1_HlP b2_HlP b3_HlP b4_HlP
b1_HlP=5.842e-16;
b2_HlP=-8.466e-09;
b3_HlP=0.06278;
b4_HlP=-4.315e+05;

%
% 1000*(b1_hlsat + b2_hlsat*(1-Trl)^(1/3) + b3_hlsat*(1-Trl)^(2/3) + b4_hlsat*(1-Trl) + b5_hlsat*(1-Trl)^(4/3));
% global b1_hgsat b2_hgsat b3_hgsat b4_hgsat b5_hgsat
b1_hgsat=-200;
b2_hgsat=+440.055;
b3_hgsat=-459.701;
b4_hgsat=+434.081;
b5_hgsat=-485.338;

%  (b1_HgP)*Pc^8 + (b2_HgP)*Pc^7 + (b3_HgP)*Pc^6 + (b4_HgP)*Pc^5 + ...
% (b5_HgP)*Pc^4 + (b6_HgP)*Pc^3 + (b7_HgP)*Pc^2 + (b8_HgP)*Pc + (b9_HgP);  
% global b1_HgP b2_HgP b3_HgP b4_HgP b5_HgP b6_HgP b7_HgP b8_HgP b9_HgP
b1_HgP=-1.223e-47; 
b2_HgP=3.658e-40; 
b3_HgP=-4.626e-33; 
b4_HgP=3.221e-26; 
b5_HgP=-1.346e-19;
b6_HgP=3.452e-13;
b7_HgP=-5.306e-07; 
b8_HgP=0.4495;
b9_HgP=-2.311e+05;

% b1_Cpg*exp(b2_Cpg*TOx) +  (b3_Cpg)*exp(b4_Cpg *TOx);
% global b1_Cpg b2_Cpg b3_Cpg b4_Cpg
b1_Cpg=308.4;
b2_Cpg=0.005439; 
b3_Cpg=6.942e-13; 
b4_Cpg=0.121;

% b1_Cpl*exp( b2_Cpl *TOx) + (b3_Cpl)*exp(b4_Cpl  *TOx);                                            % Sustancia
% global b1_Cpl b2_Cpl b3_Cpl b4_Cpl
b1_Cpl=665;
b2_Cpl=0.004254;
b3_Cpl=8.183e-10;
b4_Cpl=0.09463;

% RolOx_c*exp( b1_RolOx*(1-Trl)^(1/3) + b2_RolOx*(1-Trl)^(2/3) + b3_RolOx*(1-Trl) + b4_RolOx*(1-Trl)^(4/3) );
% global b1_RolOx b2_RolOx b3_RolOx b4_RolOx
b1_RolOx=+1.72328;
b2_RolOx=-0.83950;
b3_RolOx=+0.51060;
b4_RolOx=-0.10412;

% b1_VlP*exp((b2_VlP)*PNos1) +  (b3_VlP)*exp((b4_VlP)*PNos1);
% global b1_VlP b2_VlP b3_VlP b4_VlP b1_VgP b2_VgP b3_VgP b4_VgP
b1_VlP=0.000858;
b2_VlP=8.077e-08;
b3_VlP=9.149e-12;
b4_VlP=2.5e-06;
b1_VgP=0.09481;
b2_VgP=-1.873e-06; 
b3_VgP=0.03302;
b4_VgP=-3.344e-07;

% (b1_SgP)*Pc^7 + (b2_SgP)*Pc^6 + (b3_SgP)*Pc^5 + (b4_SgP)*Pc^4 + ...
% (b5_SgP)*Pc^3 + (b6_SgP)*Pc^2 + (b7_SgP )*Pc + (b8_SgP);
% global b1_SgP b2_SgP b3_SgP b4_SgP b5_SgP b6_SgP b7_SgP b8_SgP
b1_SgP=-1.296e-43; 
b2_SgP=3.411e-36; 
b3_SgP=-3.702e-29;
b4_SgP=2.141e-22;
b5_SgP=-7.116e-16; 
b6_SgP=1.366e-09;
b7_SgP=-0.001477;
b8_SgP=1919;

% ( b1_SlP)*Pc^3 + (b2_SlP)*Pc^2 + (b3_SlP)*Pc + (b4_SlP);
% global b1_SlP b2_SlP b3_SlP b4_SlP
b1_SlP=2.44e-18;
b2_SlP=-3.679e-11; 
b3_SlP=0.0002576;
b4_SlP=-248.2;

a_VW_NO2=5.284*10^-1;  % Atraccion entre moelesuclas de Van der Walls
b_VW_NO2=44.24*10^-6;  % Tamaño de las Moleculas Van der Walls

a_VW_N2=0.139;          %Paramtros Van der Waals Nitrogeno
b_VW_N2=0.0000391;     %Paramtros van der Waals Nitrogeno 

Tatm=273+25;

Patm=77*10^3;
g=9.8;
mmv=0.044;
mmg=0.028;
R=8.314;
RolFu=790;              % Desnidad del Etanol
Gammag=1.4;
Rv=R/mmv;
Rg=R/mmg;
RolOx_c=452;

% .---------------------------------------------------------------------
L_char=1.69;                        %%Sustancia

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
Cvg=740;
Cvw=500;
CvFu=118;

Muv=15*10^-6;
Mug=18*10^-6;
Mul_Ox=  0.000077;
Mul_Fu=  0.00110;

KlOx=102*10^-3;
Kv=19*10^-3;
Kg=26*10^-3;
KlFu=167*10^-3;

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
