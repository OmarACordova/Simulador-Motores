close all
clear all
ConstantesFluidos
R=8.314
Pv=[1*10^5:5*10^5:7*10^6]

Tvap=b1_TPvap*exp(b2_TPvap*Pv) +  b3_TPvap*exp(b4_TPvap*Pv);
Hl = (b1_HlP)*Pv.^3 + (b2_HlP)*Pv.^2 + (b3_HlP).*Pv + b4_HlP;                     % Sustancia
Hg = (b1_HgP)*Pv.^8 + (b2_HgP)*Pv.^7 + (b3_HgP)*Pv.^6 + (b4_HgP)*Pv.^5 + (b5_HgP)*Pv.^4 + (b6_HgP)*Pv.^3 + (b7_HgP)*Pv.^2 + (b8_HgP)*Pv + (b9_HgP); 
Vl = b1_VlP*exp((b2_VlP)*Pv) +  (b3_VlP)*exp((b4_VlP)*Pv);
Vg =  b1_VgP*exp((b2_VgP )*Pv) +  (b3_VgP)*exp((b4_VgP)*Pv);

figure (1)
subplot(3,1,1)
plot(Pv,Tvap )
subplot(3,1,2)
plot(Pv,Hl, Pv, Hg )
legend('Hl','Hg')
subplot(3,1,3)
plot(Pv,Vl, Pv, Vg )
legend('Vl','Vg')

Te1=[250:5:304]
Tr=Te1/304
Pvap=b1_PvapT*Te1.^3 + b2_PvapT*Te1.^2 + b3_PvapT*Te1 + b4_PvapT;
ulsat1=1000*(b1_hlsat + b2_hlsat*(1-Tr).^(1/3) + b3_hlsat*(1-Tr).^(2/3) + b4_hlsat*(1-Tr) + b5_hlsat*(1-Tr).^(4/3));
ugsat1=1000*(b1_hgsat + b2_hgsat*(1-Tr).^(1/3) + b3_hgsat*(1-Tr).^(2/3) + b4_hgsat*(1-Tr) + b5_hgsat*(1-Tr).^(4/3));   % cambiar
Cpg0 = b1_Cpg*exp(b2_Cpg*Te1) +  (b3_Cpg)*exp(b4_Cpg *Te1);
Cpl0 = b1_Cpl*exp( b2_Cpl *Te1) + (b3_Cpl)*exp(b4_Cpl  *Te1);
RolOx=RolOx_c*exp( b1_RolOx*(1-Tr).^(1/3) + b2_RolOx*(1-Tr).^(2/3) + b3_RolOx*(1-Tr) + b4_RolOx*(1-Tr).^(4/3) );

i=1;
for OFR=[1:.1:9]
% Temperatura de los Gases de combustion
    To_NOS_ETH(i) = (To_1*OFR^2 + To_2*OFR + To_3)/(OFR^3 + To_4*OFR^2 + To_5*OFR + To_6);
    % Masa Molar
    Mm_NOS_ETH(i) = (Mm_1*OFR^2 +    Mm_2.*OFR + Mm_3)/ (OFR^2 + Mm_4.*OFR + Mm_5);
    % Gamma
    Gamma_NOS_ETH(i)=(Gamma_1.*OFR^3 + Gamma_2.*OFR^2 + Gamma_3.*OFR + Gamma_4)/(OFR^2 + Gamma_5.*OFR + Gamma_6);
    %Char_Velocity(i)=(R*To_NOS_ETH(i)/Mm_NOS_ETH(i))^.5/((Gamma_NOS_ETH(i)^.5)*((2/(Gamma_NOS_ETH(i)+1))^((Gamma_NOS_ETH(i)+1)/(2*Gamma_NOS_ETH(i)-2)) ));
    BBB(i)=((Gamma_NOS_ETH(i)*Mm_NOS_ETH(i)/(R*To_NOS_ETH(i)))^.5)*((2/(Gamma_NOS_ETH(i)+1))^((Gamma_NOS_ETH(i)+1)/(2*Gamma_NOS_ETH(i)-2)));
    i=i+1
end

OFR=[1:.1:9]
figure(2)
subplot(4,1,1)
plot(Te1,Pvap )
subplot(4,1,2)
plot(Te1,ulsat1, Te1, ugsat1 )
legend('Hl','Hg')
subplot(4,1,3)
plot(Te1,Cpg0, Te1, Cpl0 )
legend('Cpg','Cpl')
subplot(4,1,4)
plot(Te1,RolOx)
legend('RolOx')


figure(3)
subplot(4,1,1)
plot(OFR,To_NOS_ETH )
legend('To')

subplot(4,1,2)
plot(OFR,Mm_NOS_ETH)
legend('Mm')

subplot(4,1,3)
plot(OFR,Gamma_NOS_ETH)
legend('Gamma')

subplot(4,1,4)
plot(OFR,BBB)
legend('BBB')