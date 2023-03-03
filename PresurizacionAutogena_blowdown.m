clear all
syms dQgwOx dQgfOx dQlfOx dQlwOx dQwatmOx 
syms dmv dTgOx dTlOx dTwOx dVgOx
syms dPv dWOx dPgOx
syms mv mgOx mlOx mwOx TgOx TlOx TgOx TwOx VgOx
syms dmLineOx VLineOx

% Varibles Definidas
syms Tf Tatm
syms PtankOx VtankOx Pv PgOx
syms hgwOx hgfOx hlwOx hlfOx hwatm
syms AgwOx Agf AlwOx Alf AwatmOx
syms Rv Rg Cvv Cvl Cvg Cvw hvap Rol mmv mmg
syms dmgOx TgpOx

% Variables de Estado
% Tg mv Vg
% Tl Tw

% conservacion de la energia
eqedo1= dQgwOx + dQgfOx - dWOx + dmv*Rv*TgOx + dmgOx*Rg*TgpOx == (mv*Cvv + mgOx*Cvg)*dTgOx % dQgw dQgf dW dmv dTg
eqedo2= dQlfOx - dQgfOx == dmv*( (Cvv+Rv)*(TgOx-Tf) - Cvl*(TlOx-Tf) + hvap )                % dQlf dQgf dmv
eqedo3= dQlwOx - dQlfOx + dWOx == dmLineOx*(0.5*VLineOx^2) + mlOx*Cvl*dTlOx                   % dQlw dQlf dW dmex? Vex dTl
eqedo4= dQwatmOx - dQlwOx - dQgwOx == mwOx*Cvw*dTwOx                                            % dQwatm dQlf dQgw dTw

dQgwOx=hgwOx*AgwOx*(TgOx-TwOx)
dQgfOx=hgfOx*Agf*(Tf-TgOx)
dQlfOx=hlfOx*Alf*(TlOx-Tf)
dQlwOx=hlwOx*AlwOx*(TwOx-TlOx)
dQwatmOx=hwatm*AwatmOx*(Tatm-TwOx)

dWOx=PtankOx*dVgOx                                                       % dVg
eqedo5= dVgOx==( ((dmv/mmv+dmgOx/mmg)*TgOx + (mv+mgOx)*dTgOx)*(Pv+PgOx) - (mgOx+mv)*TgOx*(dPv+dPgOx) )*(Rv/(Pv+PgOx)^2)                     % dmv dTg dPv

%(dQgwOx + dQgfOx - dWOx + dmv*(Rv+Cvv)*(TgOx-TgOx) + dmgOx*(Rg+Cvg)*(TgpOx-TgOx))/(mv*(Rv+Cvv) + mgOx*(Rg+Cvg)) = dTgOx % dQgw dQgf dW dmv dTg

Pex=PtankOx-Rol*0.5*VLineOx^2

% Ecuaciones de Estado
PtankOx=Pv+PgOx
Pv=Rv*TgOx*mv/VgOx
PgOx=Rg*TgOx*mgOx/VgOx
VlOx=mlOx/Rol

% Conservacion del Volumen
Vtank=VlOx+VgOx


eqedo6= dQgwFu - dQlgFu - dWFu + dmgFu*Rg*(TgpFu-TgFu) == mgFu*Cvg*dTgFu                           % dQgw dQgf dW dmv dTg
eqedo8= dQlwFu + dQlfFu + dWFu -dmLineFu*(0.5*VLineFu^2) == mlFu*CvFu*dTlFu                   % dQlw dQlf dW dmex? Vex dTl
eqedo9= dQwatmOx - dQlwOx - dQgwOx == mwOx*Cvw*dTwOx                                            % dQwatm dQlf dQgw dTw

%% Resolver

EQ=[eqedo1, eqedo2,eqedo3, eqedo4, eqedo5]
X=[dmv, dTgOx, dTlOx, dTwOx, dVgOx]

SOL=solve(EQ, X)
dmv=collect(SOL.dmv, (Cvv+Rv)*(TgOx-Tf) - Cvl*(TlOx-Tf) + hvap)
dTg=SOL.dTgOx
dTl=SOL.dTlOx
dTw=SOL.dTwOx
dVg=SOL.dVgOx

%%

dTgOx=(dQgfOx*hvap + dQgwOx*hvap - dWOx*hvap + Cvl*Tf*dQgfOx + Cvl*Tf*dQgwOx - Cvv*Tf*dQgfOx - Cvv*Tf*dQgwOx + Cvv*TgOx*dQgfOx + Cvv*TgOx*dQgwOx - Cvl*Tf*dWOx + Cvv*Tf*dWOx - Cvl*TlOx*dQgfOx - Cvl*TlOx*dQgwOx - Cvv*TgOx*dWOx + Cvl*TlOx*dWOx - Rv*Tf*dQgfOx - Rv*Tf*dQgwOx + Rv*TgOx*dQgwOx + Rv*TgOx*dQlfOx + Rv*Tf*dWOx - Rv*TgOx*dWOx + Rg*TgpOx*dmgOx*hvap + Cvl*Rg*Tf*TgpOx*dmgOx - Cvv*Rg*Tf*TgpOx*dmgOx + Cvv*Rg*TgOx*TgpOx*dmgOx - Cvl*Rg*TgpOx*TlOx*dmgOx - Rg*Rv*Tf*TgpOx*dmgOx + Rg*Rv*TgOx*TgpOx*dmgOx)/((Cvg*mgOx + Cvv*mv)*((Cvv + Rv)*TgOx + hvap - Tf*(Cvv - Cvl + Rv) - Cvl*TlOx))
%(Cvl - Cvv - Rv)*Agf*hgfOx*Tf^2 + (AgwOx*(Cvv+Rv)*hgwOx - Agf*Cvv*hgfOx )*TgOx^2

collect( (AgwOx*hgwOx*hvap - Cvv*dVgOx*(PgOx + Pv) - Rv*dVgOx*(PgOx + Pv) - Agf*hgfOx*hvap - Tf*(Agf*Cvl*hgfOx - 2*Agf*Cvv*hgfOx - AgwOx*Cvl*hgwOx + AgwOx*Cvv*hgwOx - Agf*Rv*hgfOx + AgwOx*Rv*hgwOx + Alf*Rv*hlfOx) + Agf*Cvl*TlOx*hgfOx - AgwOx*Cvl*TlOx*hgwOx - AgwOx*Cvv*TwOx*hgwOx + Cvv*Rg*TgpOx*dmgOx + Alf*Rv*TlOx*hlfOx - AgwOx*Rv*TwOx*hgwOx + Rg*Rv*TgpOx*dmgOx)*TgOx + Tf*(Cvv*dVgOx*(PgOx + Pv) - Cvl*dVgOx*(PgOx + Pv) + Rv*dVgOx*(PgOx + Pv) + Agf*hgfOx*hvap - Agf*Cvl*TlOx*hgfOx - AgwOx*Cvl*TwOx*hgwOx + AgwOx*Cvv*TwOx*hgwOx + Cvl*Rg*TgpOx*dmgOx - Cvv*Rg*TgpOx*dmgOx + AgwOx*Rv*TwOx*hgwOx - Rg*Rv*TgpOx*dmgOx) - dVgOx*hvap*(PgOx + Pv) + Cvl*TlOx*dVgOx*(PgOx + Pv) - AgwOx*TwOx*hgwOx*hvap + Rg*TgpOx*dmgOx*hvap - Cvl*Rg*TgpOx*TlOx*dmgOx + AgwOx*Cvl*TlOx*TwOx*hgwOx...
    , TgOx)

%% Presion en la camara y flujo en injectores
clear all
syms Pc 

AiOx=3.14*(0.0015*.5)^2
AiFu=3.14*(0.0015*.5)^2
NiOx=12
NiFu=3
RolOx2=900
RolFu2=900
PtankOx=4*10^6
PtankFu=4*10^6

At=3.14*(0.02*.5)^2
Gammac=1.4
Rc=8.314/0.04
Tc=3000

Vchar=(Rc*Tc*2*Gammac/(Gammac+1))^.5
RoRatio=(.5*Gammac+0.5)^(1/(1-Gammac))

eqChamber = AiOx*NiOx*(2*RolOx2*(PtankOx-Pc))^.5 + AiFu*NiFu*(2*RolFu2*(PtankFu-Pc))^.5 == At*Vchar*RoRatio*Pc/(Tc*Rc)

Pc=vpasolve(eqChamber, Pc)

mfOx=NiOx*AiOx*(2*RolOx2*(PtankOx-Pc))^.5
mfFu=NiFu*AiFu*(2*RolFu2*(PtankFu-Pc))^.5




