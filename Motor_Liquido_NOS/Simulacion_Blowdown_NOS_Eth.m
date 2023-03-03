%% Hacer Simulacion
clear all
close all
% Correr Variables fisicas generales
% MODOS DE OPERACION
% MODO 1 -VENTILACION: Valvula de ventilacion abierta apra evitar que la
% presion del tanque de Oxidante incremente
% MODO 2 -PRESURIZACION AUTOGENA: Se cierra la valvula de ventilacion para
% permitir que crezca la presion de vapor 
% MODO 3 -PRESURIZACION SOBRECARGA 
% MODO 4 -FLUJO Regulated
% MODO 5 -Flujo Blowdown

% Simulacion Blowdown: Solo con la carga inicial de Nitrogeno
% Simulacion Pressure Regula ted

ConstantesFluidos

MODO=1;
i=1;                    % Contador de tiempo de simulacion
ITotal=0;               % Impulso total acumulaod en la simulacion de la curva de motor
M_Ox=0;                 % Masa total de Oxidante expulsado
M_Fu=0;                 % Masa Total de Combustible expulsado
Dtank=6*0.0254          % Diamtro del tanque en [m] Variable de Diseño
Alg=(Dtank^2)*0.25*3.14 % Area de la interfas Liquido-Gas
Pv(1)=0.1*10^6          % Presion de los vapores de oxido nitroso Depende de la operacion
TlOx(1)=270+25;         % Temperatura de Oxidante liquido % Variable de Diseño
PTankOxReg=5*10^6       % Presion a la que el tanque de nitrogeno suministra Variable de Diseño
PTankFuReg=6*10^6       % Presion a la que el tanque de nitrogeno suministra Variable de Diseño
PTankOxRelief=6.5*10^6
mlFu(1)=2.6             % Masa de combustible % Variable de Diseño
VlOx(1)=15*10^-3        % Volumen inicial del liquido del oxidante % Variable de Diseño
UllageOx=0.5           % Variable de Diseño
UllageFu=1.1            % Variable de Diseño

MaxOpPressProp=7*10^6
MaxOpPressPres=31*10^6    % Variable de Diseño
MaxOpPressChamber=4*10^6
PgPr(1)=MaxOpPressPres         % Presion inicial del presurizante en el tanque del presuriaznante operacion
PgOx(1)=Patm            % Presion inicial del prezurizante en el tanque del Oxidante operacion
PgFu(1)=Patm            % Presion inicial del presurizante en el tanque del combustible operacion
VlFu(1)=mlFu(1)/RolFu   % Volumen inicial del liquido del combustible

% Parametros de la cara de Inyectores
LInjOx=10*10^-3;
DInjOx=2*10^-3 % Variable de Diseño
NInjOx=12; % Variable de Diseño
AInjOx=NInjOx*3.14*.25*DInjOx^2;

LInjFu=10*10^-3;
DInjFu=2*10^-3 % Variable de Diseño
NInjFu=2; % Variable de Diseño
AInjFu=NInjFu*3.14*.25*DInjFu^2

Cd=0.66

% Parametros de la Camara de combustion y Tobera
Dthroat =1*.0254;         % Variable de Diseño
DChamber=2*0.0254;          % Variable de Diseño
Athroat =3.14*0.25*Dthroat^2;
AChamber=3.14*0.25*DChamber^2;
Aexit=Athroat*5         % Variable de Diseño
Exp_Ratio=Aexit/Athroat;
L_char=1.69;
VolChamber=L_char*Athroat
LChamber=VolChamber/AChamber

% Diametro de la tuberia de la salida de lso liquidos de lso tanques
% Alimatacion del motor Operacion
WLineOx=0.889*10^-3;
DLineOx=0.0254*1/2-2*WLineOx;     % Variable de Diseño
LLineOx=1;
ALineOxo=DLineOx^2*0.25*3.14;
kLineOx=3*10^-6/DLineOx; % Rugosidad relativa de la tuberia de cobre o acero inoxidable flexible "Cooper tubing"

WLineFu=0.762*10^-3;
DLineFu=0.0254*3/8-2*WLineFu;     % Variable de Diseño
LLineFu=1;
ALineFuo=DLineFu^2*0.25*3.14
kLineFu=3*10^-6/DLineFu; % Rugosidad relativa de la tuberia de cobre o acero inoxidable flexible "Cooper tubing"

% Diametro del agujero de ventilacion del tanque, esto para permitir la
% ventilacion del tanque Operacion
DVent=0.25*0.0254;         % Variable de Diseño
AVento=0.25*3.14*(DVent)^2;

% Diametro del orificio del regualdor de presio
%Presurizacion Extra de lso tanques Operacion
DReg=(1/8)*0.0254;         % Variable de Diseño
ARego=0.25*3.14*(DReg)^2;
Cd_Press_Reg=0.5;

% Diametro del orificio de Purga
%Presurizacion Extra de lso tanques Operacion
DRelief=(1/4)*0.0254;  ;         % Variable de Diseño
ARelief=0.25*3.14*(DRelief)^2;
Cd_Relief=0.8;

% Iniciar Variables de Estado
VgOx(1)=VlOx(1)*UllageOx         
VgPr=15*10^-3;           % Variable de Diseño
VgFu(1)=VlFu*UllageFu        
TgOx(1)=Tatm
TgPr(1)=Tatm
TgFu(1)=Tatm

TlFu(1)=Tatm;

TwOx(1)=(Tatm+TlOx(1))/2;
TwFu(1)=(Tatm+TlFu(1))/2;
TwPr(1)=(Tatm+TgPr(1))/2;

TwoOx(1)=TwOx(1);
TwiOx(1)=TwOx(1);
TwoFu(1)=TwFu(1);
TwiFu(1)=TwFu(1);
TwoPr(1)=TwPr(1);
TwiPr(1)=TwPr(1);

%Variables completementarias
Trg=TgOx(1)/(273+36);            % Temepratrua reducida del los vapores del oxidnate
Trl=TlOx(1)/(273+36);            % Temepratura reducuida de los liquidos del oxidnante

% Densidad liquida del Oxidante en funcion de la temperatura
RolOx=RolOx_c*exp(b1_RolOx*(1-Trl)^(1/3) + b2_RolOx*(1-Trl)^(2/3) + b3_RolOx*(1-Trl) + b4_RolOx*(1-Trl)^(4/3)) 

%Rov(1)=(0.09481*exp((-1.873e-06 )*PNos1) +  (0.03302)*exp((-3.344e-07)*PNos1);

mv(1)=VgOx(1)*Pv(1)/(Rv*TgOx(1))            % Masa de los vapores del oxidnate
mgOx(1)=PgOx(1).*VgOx(1)/(TgOx(1)*Rg)       % Masa del presurizante en el tanque del Oxidante
mgPr(1)=PgPr(1)*VgPr/(TgPr(1)*Rg)           % Masa del presuirzante en el tanque del Presuizrante
mgFu(1)=PgFu(1)*VgFu(1)/(TgFu(1)*Rg)        % Masa del presurizante en el tanque del combustible
mlOx(1)=VlOx(1)*RolOx           % Masa del liquido de Oxidante

VtankOx=VlOx(1)+VgOx(1)         % Volumen del tanque del Oxidante
VtankFu=VlFu(1)+VgFu(1)         % Volumen del tanque del combustible

dPv=0;          %Derivada de la presion de los oxidnates tanque del oxidante
dPgOx=0;        % Derivnada de los presion del presuriante en el tanque del oxidante
PtankOx(1)=Pv(1)+PgOx(1)       % Presion total del tanque dle oxidante

tTankProp= (Dtank/2)*(3*MaxOpPressProp)/CFRP_Sigma
tTankPress=(Dtank/2)*(3*MaxOpPressPres)/CFRP_Sigma
tChamber=  (DChamber/2)*(3*MaxOpPressChamber)/CFRP_Sigma
LtankOx=(VlOx(i)+VgOx(i))/Alg
LtankFu=(VlFu(i)+VgFu(i))/Alg
LtankPr=VgPr/Alg
% Calculo de areas, coeficiente sde conveccion, largos, entalpias, etc
CalculosPropiedadesTermo

mwOx=AwatmOx*tTankProp*CFRP_Ro              % Masa de las paredes del tanque del Oxidante
mwFu=AwatmFu*tTankProp*CFRP_Ro              % Masa de las paredes del tanque del Combustible
mwPr=AwatmPr*tTankPress*CFRP_Ro             % Masa de las paredes del tanque del Presurizante
mwCh=DChamber*3.14*(LChamber+2*.25*DChamber)*tChamber*CFRP_Ro
mplumbing=2;                                    % Masa de las valvulas de llenado, alimentacion, alivio, sensores, tuberias, actuadores
M_Dry_Propulsion=mwOx+mwFu+mwCh+mplumbing
% inciializar varibles complemtarias
dQgwOx(1)=0;
dQgfOx(1)=0;
dQlfOx(1)=0;
dQlwOx(1)=0;
dQwatmOx(1)=0;

dQgwFu(1)=0;
dQlfFu(1)=0;
dQlwFu(1)=0;
dQwatmFu(1)=0;

dQgwPr(1)=0;
dQwatmPr(1)=0;

dWOx=0;
dWFu=0;
dVgFu=0;
dPOxPipe=0;
dPFuPipe=0;
dmLineOx(1)=0;                                       % Flujo de masa en la tuberias del oxidnate
dmLineFu(1)=0;                                       % Flujo de masa en la tuberias del oxidnate
VLineOx=0;
VLineFu=0;
dmgFu(1)=0;

Cvl=b1_Cpl*exp( b2_Cpl*TlOx(1)) + (b3_Cpl)*exp(b4_Cpl*TlOx(1));


dt=0.001;
time(1)=0;
dmgOx(1)=0; 
TgpOx=Tatm;

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
tiempo=0;

%%
MODO=1
while(tiempo<3*60+dt & mlOx(i)>=0 && mlFu(i)>=0)
   
    if MODO==1
        AVent=AVento;
        AReg=0;
        ALineOx=0;
        ALineFu=0;
        
        if tiempo<=0.1
            dt=0.001;
        end
        if .1<tiempo && tiempo<=3
            dt=0.01;
        end
        if 10<tiempo
            dt=.03;
        end
    end

    if MODO==2
        AVent=0;
        AReg=0;
        ALineOx=0;
        ALineFu=0;
        
        if tiempo<=0.1
            dt=0.001;
        end
        if .1<tiempo && tiempo<=3
            dt=0.01;
        end
        if 1<tiempo
            dt=.03;
        end
    end

    if MODO==3
        AVent=0;
        AReg=ARego;
        ALineOx=0;
        ALineFu=0;
        
        if tiempo<=0.1
            dt=0.001;
        end
        if .1<tiempo && tiempo<=3
            dt=0.01;
        end
        if 1<tiempo
            dt=.03;
        end
    end

    if MODO==4
        AVent=0;
        AReg=ARego;
        ALineOx=ALineOxo;
        ALineFu=ALineFuo;
        dt=0.001;
    end

    if MODO==5
        AVent=0;
        AReg=0;
        ALineOx=ALineOxo;
        ALineFu=ALineFuo;
        dt=0.003;
    end

    CalculosPropiedadesTermo

    % OXIDANTE

    dQgfOx(i)=hgfOx(i)*AgfOx*(Tf(i)-TgOx(i));  % Calor que el gas gana de la interfas
    dQlfOx(i)=hlfOx(i)*AlfOx*(TlOx(i)-Tf(i)) ; % Calor que el liquido pierde de la Interfas
    dQlwOx(i)=hlwOx(i)*AlwOx*(TwiOx(i)-TlOx(i)); % Calor que el liqudo gana de la pared
    dQgwOx(i)=hgwOx(i)*AgwOx*(TwiOx(i)-TgOx(i)); % Calor que gas Gana de la pared
    dQwatmOx(i)=hwatm*AwatmOx*(Tatm-TwoOx(i)); % Calor que la pared gana de la atmosfera
    
    PtankOx(i+1)=Pv(i)+PgOx(i);
    Pv(i+1)=Rv*TgOx(i).*mv(i)./VgOx(i);
    VlOx(i)=mlOx(i)/RolOx;   
    PgOx(i+1)=Rg*TgOx(i).*mgOx(i)./VgOx(i);  

    % Ecuaciones de Estado
    dmv = -(dQgfOx(i) - dQlfOx(i))/(hvap  + (Rv+Cvv)*(TgOx(i)-Tf(i)) + Cvl*(Tf(i)-TlOx(i)) );
    dTgOx = (dQgwOx(i) + dQgfOx(i) - dWOx + dmgOx(i)*(Rg+Cvg)*(TgpOx-TgOx(i)))/(mv(i)*(Rv+Cvv) + mgOx(i)*(Rg+Cvg));
    dTlOx = -(dmLineOx(i)*0.5*VLineOx^2 + 2*dQlfOx(i) - 2*dQlwOx(i) - 2*dWOx)/(2*Cvl*mlOx(i));
    dTwOx = -(dQgwOx(i) + dQlwOx(i) - dQwatmOx(i))/(Cvw*mwOx);
    dmVent=AVent*(PtankOx(i)/(Rv*TgOx(i)))*real(subplus((PtankOx(i)-Patm)^.5));

    if i<=1
        dVgOx = 0;
    else
        dVgOx = -(VlOx(i)-VlOx(i-1))/dt;
    end

    if  (PtankOx(i)>PTankOxRelief)
        dmOxRefiel=0.5*Cd_Relief*ARelief*(PtankOx(i)/(Rv*TgOx(i)))*real(subplus((PtankOx(i)-PTankOxReg)^.5))+ 0.5*dmOxRefiel;
    else
        dmOxRefiel=0;
    end

    %Actualizar Variables de Estado
    mv(i+1)=mv(i) + dmv*dt - dmVent*dt - dmOxRefiel*(Pv(i)/PtankOx(i))*dt;
    mlOx(i+1)=mlOx(i) - dmv*dt - dmLineOx(i)*dt;
    TgOx(i+1)=TgOx(i)+dTgOx*dt;
    TlOx(i+1)=TlOx(i)+dTlOx*dt;
    TwOx(i+1)=(TwiOx(i)+TwoOx(i))/2;
    VlOx(i+1)=VlOx(i)-dVgOx*dt;
    VgOx(i+1)=VtankOx-VlOx(i+1);
    mgOx(i+1)=mgOx(i) + dmgOx(i)*dt - dmOxRefiel*(PgOx(i)/PtankOx(i))*dt;

    TwoOx(i+1)=(Rtank*hlwOx(i)*TlOx(i)/(hlwOx(i)+Rtank)+hwatm*Tatm)./(hwatm+Rtank-Rtank^2/(hlwOx(i)+Rtank))  + dTwOx*dt;
    TwiOx(i+1)=(Rtank*hwatm*Tatm/(hwatm+Rtank)+hlwOx(i)*TlOx(i))/(hlwOx(i)+Rtank-Rtank^2/(hwatm+Rtank)) + dTwOx*dt;

    dPv=.0*(Pv(i+1)-Pv(i))/dt;
    dPgOx=.0*(PgOx(i+1)-PgOx(i))/dt;
    dWOx=PtankOx(i)*dVgOx;
    
    if mlOx(i)<=0
        break;
    end
    if mlFu(i)<=0
        break;
    end


    % PRESURIZANTE 

    dQwatmPr(i)=AwatmPr*hwatm*(Tatm-TwoPr(i));
    dQgwPr(i)=AwatmPr*hgwPr*(TwiPr(i)-TgPr(i));

    % Variables de Estado
    NgPg=mgPr(i)/mmg;
    PgPr(i+1)=NgPg*R*TgPr(i)./(VgPr-NgPg*b_VW_N2)-a_VW_N2*(NgPg/VgPr).^2;
  

    [dmgOx(i+1), TgpOx]=Pressure_Reg(PgPr(i),PTankOxReg, PtankOx(i), TgPr(i),Gammag,Rg);

    dmgOx(i+1)=0.5*dmgOx(i+1)*Cd_Press_Reg*AReg + 0.5*dmgOx(i);
    dTgPr=(dQgwPr(i)-dmgOx(i+1)*(Cvg+Rg)*TgpOx)/((Cvg+Rg)*mgPr(i));
    dTwPr=-(dQgwPr(i) - dQwatmPr(i))/(Cvw*mwPr);

    TgPr(i+1)=TgPr(i) + dTgPr*dt;
    TwPr(i+1)=TwPr(i) + dTwPr*dt;
    TwoPr(i+1)=(Rtank*hgwPr*TgPr(i)/(hgwPr+Rtank)+hwatm*Tatm)/(hwatm+Rtank-Rtank^2/(hgwPr+Rtank))  + dTwPr*dt;
    TwiPr(i+1)=(Rtank*hwatm*Tatm/(hwatm+Rtank)+hgwPr*TgPr(i))/(hgwPr+Rtank-Rtank^2/(hwatm+Rtank)) + dTwPr*dt;
    mgPr(i+1)=mgPr(i)-dmgOx(i+1)*dt-dmgFu(i)*dt;


    % COMBUSTIBLE
    dQwatmFu(i)=AwatmFu*hwatm*(Tatm-TwoFu(i));
    dQgwFu(i)=AwatmFu*hgwFu*(TwiFu(i)-TgFu(i));
    dQlwFu(i)=AlwFu*hlwFu*(TwiFu(i)-TlFu(i));
    dQlfFu(i)=Alg*hglFu(i)*(TlFu(i)-TgFu(i));

    PgFu(i+1)=mgFu(i)*Rg*TgFu(i)/VgFu(i);
    dWFu=PgFu(i)*dVgFu;

    [dmgFu(i+1), TgpFu]=Pressure_Reg(PgPr(i),PTankFuReg, PgFu(i), TgPr(i),Gammag,Rg);
    dmgFu(i+1)=0.5*dmgFu(i+1)*Cd_Press_Reg*AReg + 0.5*dmgFu(i);

    dTwFu = (-dQgwFu(i) - dQlwFu(i) + dQwatmFu(i))/(Cvw*mwFu);
    dTgFu = (dQgwFu(i) + dQlfFu(i) - dWFu + dmgFu(i+1)*Rg*(TgpFu-TgFu(i)) )/(mgFu(i)*Cvg);                     
    dTlFu =(-dmLineFu(i)*0.5*VLineFu^2 - dQlfFu(i) + dQlwFu(i)   + dWFu )/(mlFu(i)*CvFu);                   
    
    if i<=1
        dVgFu = 0;
    else
        dVgFu = -(VlFu(i)-VlFu(i-1))/dt;
    end

    if  (PgFu(i)>PTankOxRelief)
        dmFuRefiel=0.5*ARelief*(PgFu(i)/(Rg*TgFu(i)))*real(subplus((PgFu(i)-PTankFuReg)^.5)) + 0.5*dmFuRefiel;
    else
        dmFuRefiel=0;
    end

    TgFu(i+1)=TgFu(i) + dTgFu*dt;
    TlFu(i+1)=TlFu(i) + dTlFu*dt;
    TwFu(i+1)=TwFu(i) + dTwFu*dt;
    VgFu(i+1)=VgFu(i) + dVgFu*dt;

    TwoFu(i+1)=(Rtank*hgwFu*TgFu(i)/(hgwFu+Rtank)+hwatm*Tatm)/(hwatm+Rtank-Rtank^2/(hgwFu+Rtank)) + dTwFu*dt;
    TwiFu(i+1)=(Rtank*hwatm*Tatm/(hwatm+Rtank)+hgwFu*TgFu(i))/(hgwFu+Rtank-Rtank^2/(hwatm+Rtank)) + dTwFu*dt;

    mgFu(i+1)=mgFu(i) + dmgFu(i+1)*dt - dmFuRefiel*dt;
    TgFu(i+1)=TgFu(i) + dTgFu*dt;
    mlFu(i+1)=mlFu(i) - dmLineFu(i)*dt;
    VlFu(i+1)=mlFu(i)/RolFu;

    % Flujo en tuberias
    if MODO==4 || MODO==5
        Pc_guess=PtankOx(i)/2;
        Pc(i)=Pc_guess;
        for n=1:10
            Pc_guess=Pc(i);
            [dmLineOx(i+1) dmLineFu(i+1)]=MassFlow_InjLine(PtankOx(i), PgFu(i), Pc_guess, TlOx(i), Pumbling);     
            OFR(i)=dmLineOx(i+1)./dmLineFu(i+1);
              
            if OFR(i)>=9
                OFR(i)=9;
            end
            if OFR(i)<=1
                OFR(i)=1;
            end
    
            % Temperatura de los Gases de combustion
            To_NOS_ETH(i) = (-4.765e+10.*OFR(i).^2 + 1.042e+11.*OFR(i) + -1.742e+11)./(OFR(i).^3 + -1.808e+07.*OFR(i).^2 + 6.731e+07.*OFR(i) + -1.464e+08);
            % Masa Molar
            Mm_NOS_ETH(i) = (0.0292.*OFR(i).^2 + -0.06132.*OFR(i) + 0.09785)./ (OFR(i).^2 + -1.92.*OFR(i) + 4.573);
            % Gamma
            Gamma_NOS_ETH(i)=(0.004028*OFR(i).^3 + 1.185.*OFR(i).^2 + -2.961.*OFR(i) + 2.672)./(OFR(i).^2 + -2.527.*OFR(i) + 2.241);
            
            Char_Velocity(i)=(R*To_NOS_ETH(i)/Mm_NOS_ETH(i))^.5/(...
                (Gamma_NOS_ETH(i)^.5)*((2/(Gamma_NOS_ETH(i)+1))^((Gamma_NOS_ETH(i)+1)/(2*Gamma_NOS_ETH(i)-2)) ));
            BBB=((Gamma_NOS_ETH(i)*Mm_NOS_ETH(i)/(R*To_NOS_ETH(i)))^.5)*((2/(Gamma_NOS_ETH(i)+1))^((Gamma_NOS_ETH(i)+1)/(2*Gamma_NOS_ETH(i)-2)));
            Char_G(i)=Pc_guess*BBB;
            Pc_guess=(dmLineOx(i+1)+dmLineFu(i+1))/(BBB*Athroat);
            Pc(i)=Pc_guess*.2 + Pc(i)*.8;
        end
    
        VLineOx=dmLineOx(i+1)/(ALineOx*RolOx)
        VLineFu=dmLineFu(i+1)/(ALineFu*RolFu);
        
        % Para encontrar el numero de Mach y la presion de salida de tobera
        % usaremos las ecuaciones de flujo isentropico. Ya que no podemos
        % despejar el Mach de la ecuacion, pero si conocemos el expanction ratio. Propondremos Mach
        % y calcularemos el expantion ratio correspondiente a este Mach,
        % hasta que el Expnation ratio calculado sea igual al Real
        Mach_Guess=1;
        Exp_Ratio_Guess=1;
        while Exp_Ratio_Guess<Exp_Ratio
            Exp_Ratio_Guess=(1/Mach_Guess)*((2/(Gamma_NOS_ETH(i)+1))*(1+0.5*(Gamma_NOS_ETH(i)-1)*Mach_Guess^2))^((Gamma_NOS_ETH(i)+1)/(2*Gamma_NOS_ETH(i)-2));
            Mach_Guess=Mach_Guess+0.05;
        end
    
        Mach(i)=Mach_Guess;
        Pexit(i)=Pc(i)*(1+.5*(Gamma_NOS_ETH(i)-1)*Mach(i)^2)^(-Gamma_NOS_ETH(i)/(Gamma_NOS_ETH(i)-1));
        Thrust(i)=Athroat*Pc(i)*( ((2*Gamma_NOS_ETH(i).^2)/(Gamma_NOS_ETH(i)-1)) * ((2/(Gamma_NOS_ETH(i)+1))^((Gamma_NOS_ETH(i)+1)/(Gamma_NOS_ETH(i)-1)) )*...
            ( 1-(Pexit(i)/Pc(i))^( (Gamma_NOS_ETH(i)-1)/Gamma_NOS_ETH(i) ))   )^.5 + Aexit*(Pexit(i)-Patm);
        
        ITotal=ITotal+Thrust(i)*dt;
        M_Ox=M_Ox+dmLineOx(i)*dt;
        M_Fu=M_Fu+dmLineFu(i)*dt;
        Isp(i)=0.1*Thrust(i)/(dmLineOx(i)+dmLineFu(i));
    else
        Thrust(i)=0;
        Isp(i)=0;
        Pc(i)=0;
        Pexit(i)=0;
        dmLineOx(i+1)=0;
        dmLineFu(i+1)=0;     
        OFR(i)=0;
    end

    time(i+1)=time(i)+dt;
    i=i+1;
    tiempo=tiempo+dt;
    [i,tiempo]

end

%%
close all

figure(1)
title('Temperatura sustncias[K] vs tiempo [s]')
hold on
plot(time(1:end-1), TgOx(1:end-1), 'LineWidth', 2.5)
plot(time(1:end-1), Tf, [0,time(end)], [Tatm, Tatm])
plot(time(1:end-1), TlOx(1:end-1), 'LineWidth', 2.5)
plot(time, TgPr, 'LineWidth', 2.5)
plot(time, TgFu, 'LineWidth', 2.5)
plot(time, TlFu, 'LineWidth', 2.5)
hold off
legend('Mix Ox','L-VapInterf-Ox','Atm','LOx', 'Pr', 'FuG', 'FuL' )

figure(2)
title('Temperatura Paredes[K] vs tiempo [s]')
hold on 
plot(time, TwOx, time, TwoOx, time, TwiOx)
plot(time, TwPr, time, TwoPr, time, TwiPr)
plot(time, TwFu, time, TwoFu, time, TwiFu)
hold off
legend('WOx','WoOx','WiOx', 'WPr','WiPr','WoPr', 'WFu','WiFu','WoFu')

figure(3)
title('Presiones [Pa] vs tiempo [s]')
hold on 
plot(time(1:end-1), Pv(1:end-1)+PgOx(1:end-1), 'LineWidth', 2) 
plot(time(1:end-1), Pv(1:end-1), 'LineWidth', 2)
plot(time(1:end-1), Pvapl, time, PgOx)
plot(time(1:end-1), PgPr(1:end-1), 'LineWidth', 2)
plot(time(1:end), PgFu, 'LineWidth', 2)
hold off
legend('Press Mix','Press VapOx', 'Press Vap Temp Liq', 'Press Vap Temp Gas', 'Press N2 Ox', 'Press Press', 'Press FuG')

figure(4)
title('Flujo de Calor [Watt] vs tiempo [s]')
hold on
plot(time(1:end-2), dQgwOx(1:end-1), time(1:end-2), dQlwOx(1:end-1), time(1:end-2), dQgfOx(1:end-1), time(1:end-2), dQlfOx(1:end-1))
plot(time(1:end-2), dQwatmOx(1:end-1))
plot(time(1:end-1), dQwatmPr, time(1:end-1), dQgwPr)
plot(time(1:end-2), dQgwFu(1:end-1), time(1:end-2), dQlwFu(1:end-1), time(1:end-2), dQlfFu(1:end-1))
plot(time(1:end-2), dQwatmFu(1:end-1))
hold off
legend('Mix-WallOx','Liq-WallOx','Mix-Inter Ox','Liq-Inter Ox','Wall-Atm Ox', ...
    'Wall-Atm Pr', 'Wall-Pr' ,...
    'N2-WallFu', 'Liq-WallFu', 'N2-Ful', 'Wall-Atm Fu ')

figure(5)
title('Masas [Kg] vs tiempo [s]')
hold on 
plot(time(1:end-1), mv(1:end-1), 'LineWidth', 2)
plot(time(1:end-1), mlOx(1:end-1), 'LineWidth', 2)
plot(time(1:end-1), mlFu(1:end-1), 'LineWidth', 2)
plot(time, mgOx, time, mgPr, time, mgFu)
hold off
legend('Mass NO2 Vap', 'Mass NO2 Liq', 'Mass Fu Liq','Mass N2 Ox', 'Mass N2 Press', 'Mass N2 Fu')

figure(6)
title('Volumen [m3] vs tiempo [s]')
hold on
plot(time, VgOx, time, VlOx)
plot(time, VgFu, time, VlFu)
hold off
legend('Vol Vap Ox', 'Vol Liq Ox','Vol Vap Fu', 'Vol Liq Fu')

figure(9)
plot(time(1:end-1),hgfOx,time(1:end-1),hlfOx,  time(1:end-1),hlwOx,time(1:end-1),hgwOx,  time(1:end-1),hglFu)
title('Coeficientes de Conveccion Ox [W/K*m2]')
legend('Vap-L Ox','L-Vap Ox','L-Wall Ox', 'Vap-Wall Ox', 'G-L Fu')

figure(7)
subplot(2,2,1)
plot(time(1:end-1), Thrust(1:end) )
title('Empuje [N]')

subplot(2,2,2)
plot(time(1:end), dmLineOx(1:end), time(1:end), dmLineFu(1:end), time(1:end), dmgOx(1:end), time(1:end), dmgFu(1:end) )
title('Flujos de Masa [kg/s]')
legend('mlOx','mlFu', 'mgPrOx', 'mgPrFu')

subplot(2,2,3)
plot(time(1:end-1), 100*OFR(1:end), time(1:end-1), Isp(1:end))
title('OFR & ISP')
legend('OFR','ISP')

subplot(2,2,4)   
plot(time(1:end-1), Pc(1:end), time(1:end), PtankOx(1:end), time(1:end), PgFu(1:end), time(1:end-1), Pexit )
title('Presion Camara Combustions [MPa]')
legend('Pc', 'PtankOx', 'PgFu', 'PExit')

figure(8)
subplot(3,1,1)
hold on
scatter( PtankOx(end-15*10^3:end)-Pc(end-15*10^3:end), mlOx(end-15*10^3:end) )
scatter( PgFu(end-15*10^3:end)-Pc(end-15*10^3:end), mlFu(end-15*10^3:end)  )
hold off
xlim([0 2*10^6])
title('Flujo de Masa vs dP')
legend('mlOx','mlFu')

subplot(3,1,2)
hold on
scatter(PtankOx(end-15*10^3:end), mlOx(end-15*10^3:end))
scatter(   PgFu(end-15*10^3:end), mlFu(end-15*10^3:end))
hold off
xlim([0 5*10^6])
title('Flujo de Masa vs Ptank')
legend('mlOx','mlFu')

subplot(3,1,3)
hold on
scatter(Pc(end-15*10^3:end), mlOx(end-15*10^3:end))
scatter(Pc(end-15*10^3:end), mlFu(end-15*10^3:end))
hold off
xlim([0 3*10^6])
title('Flujo de Masa vs Pchamber')
legend('mlOx','mlFu')

figure(10)
subplot(2,2,1)
plot(time(1:end-1), To_NOS_ETH)
title('Temperatura Combustion')

subplot(2,2,2)
plot(time(1:end-1), Mm_NOS_ETH)
title('Masa molar Combustion')

subplot(2,2,3)
plot(time(1:end-1), Gamma_NOS_ETH)
title('Gamma combustion')

subplot(2,2,4)
plot(time(1:end-1), Char_Velocity)
title('Velocidad Caracteritica ')

Impulso=ITotal(end)
M_Ox(end)
M_Fu(end)


[time'*1000, Pv'/1000, PgOx'/1000, mv', mlOx', TgOx', TlOx', TwOx', VgOx'*1000, VlOx'*1000, ];

%Regulated 
% PReg=4MPa N2=20L=15MPa
%Itotal=1.2609e+05 Temperatura=25°C Ullage=30% Etnol/NOS=7/40  

% PReg=5MPa N2=10L=15MPa
%Itotal=1.2774e+05 Temperatura=25°C Ullage=30% Etnol/NOS=7/40  

% PReg=5MPa N2=15L=25MPa
%Itotal=1.2992e+05 Temperatura=25°C Ullage=30% Etnol/NOS=7/40  PReg=4MPa

%PReg=5.5MPa N2=15L=25MPa
%Itotal=1.2715e+05 Temperatura=25°C UllageOx=20% UllageOx=100% Etnol/NOS=7/40 
% 

%PReg=5.5MPa N2=15L=25MPa
%Itotal=1.2937e+05 Temperatura=25°C UllageOx=20% UllageOx=100% Etnol/NOS=7/40 
%Con presion atmosferica 77Kpa


% Blowdown
%PReg=4MPa 
%Itotal=1.0722e+05 Temperatura=25°C Ullage=30% Etnol/NOS=7/40
%sobro mucho etanol y la ezcla fue muy rica en oxigeno, necesita una manera de presurisarse

%PReg=5MPa 
%Itotal=9.9794e+04 Temperatura=25°C UllageOx=10% UllageOx=100% Etnol/NOS=7/40 
% sobro mucho NOS

%PReg=5MPa 
%Itotal=1.2207e+05 Temperatura=25°C UllageOx=10% UllageOx=70% Etnol/NOS=7/40 
% sobro mucho NOS

%PReg=5MPa 
%Itotal=1.2114e+05 Temperatura=25°C UllageOx=50% UllageOx=100% Etnol/NOS=7/40 
% sobro etaanol

%PReg=5MPa 
%Itotal=1.2460e+05 Temperatura=25°C UllageOx=30% UllageOx=100% Etnol/NOS=7/40 
% sobro etaanol

%PReg=5.5MPa 
%Itotal=1.1923e+05 Temperatura=25°C UllageOx=20% UllageOx=100% Etnol/NOS=7/40 
% sobro poco NOS





