Hn2=py.CoolProp.CoolProp.PropsSI('H','T',273,'P',20*10^6,'N2')
Dn2=py.CoolProp.CoolProp.PropsSI('D','T',273,'P',20*10^6,'N2')
Hhe=py.CoolProp.CoolProp.PropsSI('H','T',273,'P',20*10^6,'He')
Dhe=py.CoolProp.CoolProp.PropsSI('D','T',273,'P',20*10^6,'He')
% i=1;
% for D=60:-1:1
%     d(i)=D;
%     Pn2(i)=py.CoolProp.CoolProp.PropsSI('P','D',D*8,'H',Hn2,'N2');
%     Tn2(i)=py.CoolProp.CoolProp.PropsSI('T','D',D*8,'H',Hn2,'N2');
%     Phe(i)=py.CoolProp.CoolProp.PropsSI('P','D',D,'H',Hhe,'He');
%     The(i)=py.CoolProp.CoolProp.PropsSI('T','D',D,'H',Hhe,'He');
%     i=i+1;
% end
% figure(1)
% plot3(d,Pn2,Tn2,d,Phe,The)
% legend('N2','He2')
% xlabel('Ro (kg/m3)')
% ylabel('P(Pa)')
% zlabel('T(K)')

[T_mix, Ro_1, Ro_2]=Mix_Equilibrium_1(Hhe+Hn2, 1,1,.1,'He','N2')


