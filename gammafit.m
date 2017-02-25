function [fit] = gammafit(xdata,ydata) %x is time
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
figure(4000)
plot(xdata(1:endgamma),lngamma,xdata(1:endgamma),predictln)

%fitting only from the start point on...
%finding the start point:
starttimeonwards = xdata(xdata>=0);

%now trimming away the y values that exist before 0
xvectorlength = length(starttimeonwards);
yvectorlength = length(ydata);
i = yvectorlength-xvectorlength;
newyrange = ydata(i+1:yvectorlength);

%making nice names for shit

x = starttimeonwards;
y = newyrange;
%disp(length(x))
%disp(length(y))

%fitting the data
usnume = (x.*y);
usdeno = (x.*x);

nume = sum(usnume);
deno = sum(usdeno);

anot = nume/deno;

%anot is the same thing as -t/tau, so...
predictedlnGamma = x.*anot; %this is the slope

fit = [x,predictedlnGamma];



end

