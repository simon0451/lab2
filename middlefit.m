function [output] = middlefit(xdata,ydata) %x input is time, y input is temperature from data
%the inputs are shifted such that the 0 value of the xdata is the beginning
%of the event.

%this function provides an output array containing the information to be
%plotted on the middle plot of the 3-plot subplot

gamma=(ydata(end)-ydata)/(ydata(end)-ydata(1));
for i=1:length(xdata)
    if gamma(i)<0.05
        endgamma=i;
        break
    end
end

lngamma=log(gamma(1:endgamma));
num=lngamma.*xdata(1:endgamma);
den=xdata(1:endgamma).^2;
ao=sum(num)/sum(den);


Tfinal = ydata(end); %the final temperature of the data

Tinitloc = find(xdata==0); %the location in the array where time is zero

Tinitial = ydata(Tinitloc); %the Temperature value where time is zero

tau = (1/ao)*-1;

Tpred = Tfinal - (Tfinal-Tinitial)*exp(-xdata/tau);


output = [xdata,ydata,xdata,Tpred]; %[data time, data Temperature, predicted time (same as data time), predicted Temperature]


end