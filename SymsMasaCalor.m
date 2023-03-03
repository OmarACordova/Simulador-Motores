close all
syms nga1 nga2 nla1 nla2 dngb ngb2 nlb2                     %Moles de gases y liquidos 
syms uga1 ula1 hgb Cvla Cvga Cplb Cpgb                      %Energias Internas y Entalpias 
syms Qlg hvapa hvapb                                        %Tranferencia de Calor Vapor liquido
syms Tmix Tvapa Tvapb Tl                                    %Temperaturas de Transferencia de calor Vapor-Liquido
syms Qwg Qwl Alg Agw Alw Twgi Twgo Twli Twlo Tatm Kw tw     %Tranferecia de calor con las paredes
syms hlw hgw hatmw hlg                                      %Coeficientes de conveccion
syms Ptank1 Ptank2 Pa1 Pa2 Pb2 Ta1 R Vg
%syms Pb0 Tb0 Gamma
syms dt time
syms Trefvapb Trefvapa Pref

%Variables conocidas
% nga1 nla1 dngb                            % nga2 nla2 ngb2 nlb2
% uga1 ula1 hgb Cvla Cvga Cplb Cpgb         %
% hvapa hvapb                               % Qlg
%                                           % Tmix Tvapa Tvapb Tl
% Alg Agw Alw Tatm Kw tw                    % Qwg Qwl Twgi Twgo Twli Twlo   
% hlw hgw hatmw hlg                         %
% Pa1 Ta1 R Vg                              % Ptank2 Pa2 Pb2
% dt time
% Trefvapb Trefvapa Pref
%Ecuaciones de balance de energia
nga2=((nga1*uga1+dt*dngb*hgb-Qlg+Qwg)/Tmix - ngb2*Cpgb)/Cvga            % de aqui se despeja Tmix
nla2=((nla1*ula1+Qlg+Qwl)/Tl - nlb2*Cplb)/Cvla
Qlg=(nla2-nla1)*(Cvga*(Tmix-Tvapa)+hvapa+Cvla*(Tvapa-Tl))...
   +(nlb2-nlb1)*(Cpgb*(Tmix-Tvapb)+hvapb+Cplb*(Tvapb-Tl))
Qlg=Alg*hlg*(Tmix-Tl)

Qwg=Agw*hgw*(Tmix-Twgi)
Qwg=Agw*hatmw*(Twgo-Tatm)
Qwg=(Agw*Kw/tw)*(Twgo-Twgi)

Qwl=Alw*hlw*(Tl-Twli)
Qwl=(Alw*Kw/tw)*(Twli-Twlo)
Qwl=Alw*hatmw*(Twlo-Tatm)
%Twlo=(hatmw*Tatm*tw/Kw+Twli)/(1-hatmw*tw/Kw)
%hatmw*Twlo-hatmw*Tatm=(Kw/tw)*Twli-(Kw/tw)*Twlo
%Twlo=(Twli-Tatm*hatmw*tw/Kw)/(hatmw*tw/Kw +1)

%Ecuaciones de balance de masa
nga1=nga2+nla2 - nla1
ngb2=time*dngb - nlb2

%Entradas de los tanues presurizantes a los tanques de propelente
%Tb1=Tb0*(Ptank/Pb0)^((Gamma-1)/Gamma)               %Temperatura de los gases presuruzados despues de expandirse 

%Ecuaciones de Estado, Gas ideal
nga1=Pa1*Vg/(Ta1*R)         %estado incial de gas del propelente solo en equilibrio
nga2=Pa2*Vg/(Tmix*R)        %Estado siguiente de gas del propelente

%ngb1=Pb1*Vg/(Tb1*R)
ngb2=Pb2*Vg/(Tmix*R)

Ptank1=Pa1
Ptank2=Pa2+Pb2

Tvapa=(1/Trefvapb-log(Ptank/Pref)*R/hvapa)^-1        %Clausuis clayperion para a
Tvapb=(1/Trefvapa-log(Ptank/Pref)*R/hvapb)^-1        %Clausuis clayperion para b

% Ecuaciones de conservacion del Volumen

%Ecuaciones dinamicas
Paeq=Pavap*(nla2/(nlb2+nla2))              %Ley de Raoult
Pbeq=Pbvap*(nlb2/(nlb2+nla2))              %Ley de Raoult
%%
