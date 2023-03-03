close all
ConstantesFluidos

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
