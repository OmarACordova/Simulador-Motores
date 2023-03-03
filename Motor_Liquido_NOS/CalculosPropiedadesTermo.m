% Gemetria de los Tanques


AgfOx=Alg;
AlfOx=Alg;

AgwOx=Dtank*3.14*LtankOx*(VgOx(i)/VtankOx) + Alg;
AlwOx=Dtank*3.14*LtankOx*(VlOx(i)/VtankOx) + Alg;
AwatmOx=AlwOx + AgwOx;

AgwPr=Dtank*3.14*LtankPr + 2*Alg;
AwatmPr=AgwPr;

AgwFu=Dtank*3.14*LtankFu*(VgFu(i)/VtankFu) + Alg;
AlwFu=Dtank*3.14*LtankFu*(VlFu(i)/VtankFu) + Alg;
AwatmFu=AlwFu + AgwFu;

% Propiedades de los fluidos
RolOx=RolOx_c*exp(b1_RolOx*(1-Trl)^(1/3) + b2_RolOx*(1-Trl)^(2/3) + b3_RolOx*(1-Trl) + b4_RolOx*(1-Trl)^(4/3));
ugsat1=1000*(b1_hgsat + b2_hgsat*(1-Trg)^(1/3) + b3_hgsat*(1-Trg)^(2/3) + b4_hgsat*(1-Trg) + b5_hgsat*(1-Trg)^(4/3));
ulsat1=1000*(b1_hlsat + b2_hlsat*(1-Trl)^(1/3) + b3_hlsat*(1-Trl)^(2/3) + b4_hlsat*(1-Trl) + b5_hlsat*(1-Trl)^(4/3));
hvap=ugsat1-ulsat1;

Cvl=b1_Cpl*exp( b2_Cpl *TlOx(1)) + (b3_Cpl)*exp(b4_Cpl *TlOx(1));
Cvv = b1_Cpg*exp(b2_Cpg*TgOx(i)) +  (b3_Cpg)*exp(b4_Cpg *TgOx(i));

Pvapl(i)=b1_PvapT*TlOx(i).^3 + b2_PvapT*TlOx(i).^2 + b3_PvapT*TlOx(i) + b4_PvapT;

Cvg=740;
Cvw=500;
CvFu=118;

Muv=15*10^-6;
Mug=18*10^-6;
Mul_Ox=  0.000077;

KlOx=102*10^-3;
Kv=19*10^-3;
Kg=26*10^-3;
KlFu=167*10^-3;


Ktank=15;
Rtank=tTankProp*Ktank;

Beta_lOx=1/128;
Beta_lFu=1/912;

Rov=mv(i)./VgOx(i);

Cpvg=((Cvv+Rv)*mv(i)+(Cvg+Rg).*mgOx(i))./(mgOx(i)+mv(i));
Kvg=(Kv*mv(i)+Kg.*mgOx(i))./(mgOx(i)+mv(i));
Muvg=(Muv*mv(i)*mmv+Mug.*mgOx(i)*mmg)./(mgOx(i)*mmv+mv(i)*mmg);
Rvg=(Rv*mv(i)*mmv+Rg.*mgOx(i)*mmg)./(mgOx(i)*mmv+mv(i)*mmg);

a=551.8;
b=-2.183e+05;
c=+ 2.174e+07 - Pv(i);

if Pv(i)>0.5*10^6
    Tf(i)=(-b+(b^2-4*a*c)^.5)/(a*2);
else
    Tf(i)=273-50;
end

% NOSLiq - Interfaz Numeros adimensionales
PrlOx=Cvl*Mul_Ox/KlOx;
GrlfOx=g*Beta_lOx*abs(Tf(i)-TlOx(i))*(.25*Dtank)^3/(Mul_Ox/RolOx)^2;
RalfOx=PrlOx*GrlfOx;
NulfOx=0.14*(RalfOx)^.333;
hlfOx(i)=NulfOx*KlOx/(0.25*Dtank);

% NOSLiq - Pared Numeros adimensionales
GrlwOx=g*Beta_lOx*abs(TwiOx(i)-TlOx(i))*(VlOx(i)/Alg)^3/(Mul_Ox/RolOx)^2;
RalwOx=PrlOx*GrlwOx;
NulwOx=.59*RalwOx^.25;
hlwOx(i)=NulwOx*KlOx/(VlOx(i)/Alg);

% NOSVap - Interfaz Numeros adimensionales
PrgOx=Cpvg*Muvg/Kvg;
GrgfOx=g*(1/TgOx(i))*abs(TgOx(i)-Tf(i))*(0.25*Dtank)^3/( (Muvg*Rvg*TgOx(i)./PtankOx(i))^2 );
RagOx=PrgOx*GrgfOx;
NugfOx=0.15*(RagOx)^.333;
hgfOx(i)=NugfOx*Kvg/(0.25*Dtank);

% NOSVap - Pared 
GrgwOx=g*(1/TgOx(i))*abs(TwiOx(i)-TgOx(i))*(VgOx(i)/Alg)^3/( (Muvg*Rvg*TgOx(i)./PtankOx(i))^2 );
RagwOx=PrgOx*GrgwOx;
NugwOx=.59*RagwOx^.25;
hgwOx(i)=NugwOx*Kvg/(VgOx(i)/Alg);

% Etanol - Interfas Liq
PrlFu=CvFu*Mul_Fu/KlFu;
GrlgFu=g*Beta_lFu*abs(TgFu(i)-TlFu(i))*(0.25*Dtank)^3/(Mul_Fu/RolFu)^2;
RalgFu=PrlFu*GrlgFu;
NulgFu=0.59*RalgFu^.25;

% Etanol - Interfas Gas
PrgFu=Cvg*Mug/Kg;
GrglFu=g*(1/TgFu(i))*abs(TlFu(i)-TgFu(i))*(0.25*Dtank)^3/(Mug*Rg*TgFu(i)./PgFu(i))^2;
RaglFu=PrgFu*GrglFu;
NuglFu=0.59*RaglFu^.25;

hglFu(i)=max(NulgFu*KlFu, NuglFu*Kg)/(Dtank*.25);

% Etanol - N2 Pared 
GrgwFu=g*(1/TgFu(i))*abs(TwiFu(i)-TgFu(i))*(VgFu(i)/Alg)^3/(Mug*Rg*TgFu(i)./PgFu(i))^2;
RagwFu=PrgFu*GrgwFu;
NugwFu=0.59*RagwFu^.25;
hgwFu=NugwFu*Kg/(VgFu(i)/Alg);

% Etanol - Liq Pared 
GrlwFu=g*Beta_lFu*abs(TwiFu(i)-TlFu(i))*(VlFu(i)/Alg)^3/(Mul_Fu/RolFu)^2;
RalwFu=PrlFu*GrlwFu;
NulwFu=0.59*RalwFu^.25;
hlwFu=NulwFu*KlFu/(VlFu(i)/Alg);


% N2 - ParedN2 Numeros adimensionales
PrPr=Cvg*Mug/Kg;
GrPr=g*(1/TgPr(i))*abs(TwiPr(i)-TgPr(i))*LtankPr^3/(Mug*Rg*TgPr(i)/PgPr(i))^2;
RaPr=PrPr*GrPr;
NuPr=0.59*RaPr^.25;
hgwPr=NuPr*Kg/LtankPr;

hwatm=50;