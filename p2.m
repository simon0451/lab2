function [output] = p2(xdata,ydata)
%this function outputs an array with ready to plot info for the two curves
%and the residuals

%where xdata is the time and ydata is the temperature

Tfinal = ydata(end); %the final temperature of the data
Tinitloc = find(xdata==0); %the location in the array where time is zero
Tinitial = ydata(Tinitloc); %the Temperature value where time is zero
Tattau = Tinitial+.632*(Tfinal - Tinitial); %this is the temperature at the point where t is equal to tau
%disp(Tattau)

%find the time at which the temperature is equal to Tattau - this is our
%time constant tau

%conditional statement 
if Tfinal<Tinitial
    tauindex = find(ydata<=Tattau,1); %the location of tau in the array
else
    tauindex = find(ydata>=Tattau,1);
end
    

tau = xdata(tauindex);

Tpred = Tfinal-(Tfinal-Tinitial)*exp(-xdata/tau);

residuals = ydata-Tpred;

output = [xdata,ydata,xdata,Tpred,xdata,residuals]; %[data time, data Temperature, predicted time (same as data time), predicted Temperature]

end