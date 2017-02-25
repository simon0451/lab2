function [logGamma] = timecon(Tarray,Tfinal)
%Tinitial is assumed to be the first point in the array
Tkel = Tarray;%+273.15;
Tfinal = Tfinal;%+273.15;
Gamma = (Tfinal-Tkel)./(Tfinal-Tkel(1));
disp(Gamma(1))
disp(Gamma(2000))
disp(Gamma(end))
logGamma = log(Gamma);



end

