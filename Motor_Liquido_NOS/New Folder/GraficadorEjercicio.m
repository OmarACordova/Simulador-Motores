close all
clear all
G_SPI(1)=1
j=1;
for Pt=1*10^6:5*10^5:7*10^6
    i=1;
    for Pc=1*10^6:1*10^5:7*10^6
        [G_SPI(j,i), G_HEM(j,i)]=FlujoMasa(Pt,Pc);
        i=i+1
    end
    j=j+1;
end

[Pt, Pc]=meshgrid(1*10^6:5*10^5:7*10^6, 1*10^6:1*10^5:7*10^6);
figure(1)
plot3(Pt,Pc,G_SPI)
title("G-SPI")
xlabel("Pt")
ylabel("Pc")
zlabel("G_SPI")

figure(2)
plot3(Pt,Pc,G_HEM)
title("G-HEM")
xlabel("Pt")
ylabel("Pc")
zlabel("G_HEM")

