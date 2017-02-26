function [output] = Syx(xdata,ydata)
%Computers the Syx of a data set

tstart = find(xdata==0,1);
xnew = xdata(tstart:length(xdata));
ynew = ydata(tstart:length(ydata));

fitdata = polyfit(xnew,ynew,1);
yfromline = fitdata(1)*xnew+fitdata(2);

Yc = yfromline; %The value of y predicted by the polynomial equation for a given value of x
tvp = 2.262; %For N = 10, 95% confidence
yiyci = (ynew-yfromline).^2;
sumyiyci = sum(yiyci);
Syx = (sumyiyci/(length(ynew)-1))^.5; %standard error of the fit

output = Syx;

end
