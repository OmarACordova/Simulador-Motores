% Esta funcion entrega mass Flux para una presion tanque Pt y una presion
% de camara de combustion Pc
function [G_SPI, G_HEM]=FlujoMasa(Pt,Pc)
    Propiedades
    if Pt>Pc
        G_SPI=(2*(Pt-Pc)/Vl)^.5
        G_HEM=(Pt+Pc)^.5
    else
        G_SPI=0
        G_HEM=0
    end
end
