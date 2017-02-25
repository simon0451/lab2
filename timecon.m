function [logGamma] = timecon(Tarray,Tfinal)
%this funtion produces the ln of Gamma
%the function name is a misnomer

%Tinitial is assumed to be the first point in the array
Tkel = Tarray;
Gamma = (Tfinal-Tkel)./(Tfinal-Tkel(1));
logGamma = log(Gamma);

output = [logGamma];



end

