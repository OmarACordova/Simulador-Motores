close all
clear all
%
% Este programa nos calculara el comportameinto de una sola sustencia
% compuesta por vapor y liquido que no se encuentra en equilibrio
% y que legran a equilibrio debido a la conveccion entre las interfaz
% liquido vapor
%

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

% Varibles de entrada       % Var. Faciles calcular     % Varible desconocidas 
%                           % nga1 nla1                 % nga2 nla2
% uga1 ula1 Cvla Cvga       %                           % uga2 ula2
% hvapa Alg Qlg1            %                           % Qlg2 dqlg
% Ta1 Tl1 Pa1               % Tvapa1                    % Ta2 Tvapa2 Tl2
% hlg
% Pa1 R Vg1 Vl1 mma         % Vtank                     % Pa2 Vg2 Vl2
% Rola
% dt
% Tvaparef Pref

%nga2, Ta2, dqlg, nla2, Tl2, Vl2 
%Conseervacion de al energia
eq1=nga2*(uga1+Cvla*(Ta2-Ta1))==(nga1*uga1-dqlg)            
eq2=nla2*(ula1+Cvga*(Tl2-Tl1))==(nla1*ula1+dqlg)
eq3=dqlg==dt*(nla2-nla1)*(Cvga*(Ta1-Tvapa1)+hvapa+Cvla*(Tvapa1-Tl1))
eq4=dqlg==dt*Alg*hlg*(Ta1-Tl1)
%conservacion de al masa
eq5=nga1==nga2+nla2-nla1
%Convervacion del Volumen
eq6=Vl1==Vl2+Vg2-Vg1



%%
%EQ=[eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13,eq14,eq15];
%X=[nga2,nla2,uga2, ula2, Qlg, Ta2, Tvapa2, Tl2, Pa2, Vg2, Vl2];

EQ=[eq1,eq2, eq3, eq4, eq5, eq6 ];
X=[nga2, Ta2, dqlg, nla2, Tl2, Vl2  ];

Sol=solve(EQ, X)

nga2=simplify(Sol.nga2)
nla2=simplify(Sol.nla2)
dqlg=simplify(Sol.dqlg)
Ta2=simplify(Sol.Ta2)
Tl2=simplify(Sol.Tl2)
Vl2=simplify(Sol.Vl2) 

