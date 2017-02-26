function [output] = Syx(xdata,ydata)
%Computers the Syx of a data set
fitdata = polyfit(xdata,ydata,1);
yfromline = fitdata(1)*xdata+fitdata(2);

Yc = yfromline; %The value of y predicted by the polynomial equation for a given value of x
tvp = 2.262; %For N = 10, 95% confidence
yiyci = (ydata-yfromline).^2;
sumyiyci = sum(yiyci);
Syx = (sumyiyci/(length(ydata)-1))^.5; %standard error of the fit

output = Syx;

end
