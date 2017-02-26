function [output] = gammafit(xdata,ydata) %x is time
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
predictln=ao*xdata(1:endgamma);
% figure(4000)
% plot(xdata(1:endgamma),lngamma,xdata(1:endgamma),predictln)

output = [xdata(1:endgamma),lngamma,xdata(1:endgamma),predictln];



end

