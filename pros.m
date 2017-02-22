function [output] = pros(time,involtage,startingtime)
offset = 0;

% %smooth data, locate peaks and valleys
% zerovoltage = involtage-offset; %this allows the data to be offset if desired (good for peakdet)
% span=50; %size of the averaging window
% mask=ones(span,1)/span;
% voltage=conv(zerovoltage,mask,'same'); %finds moving average and cleans up noise
voltage = smooth(involtage,51);
%cleans up sharp rise and fall spikes in data caused by the convolution
%function
% b = span/2;
% e = length(convoltage)-b;
% vout = convoltage(b:e);
% tunshifted = newtime(b:e);
% tshift = 0-tunshifted(1);
% tout = tunshifted+tshift;

for i=1:length(time)
    if time(i)>startingtime
        basetime=i;
        break
    end
end
baseline=mean(voltage(1:basetime)); %average of data in baseline region
basedev=std(voltage(1:basetime)); %standard deviation of baseline region
threshold=5*basedev; %threshold to define the start of an event - 5 sigma

for i=1:length(time)
    if (abs(voltage(i)-baseline)>threshold)
        starttime=i;
        break
    end
end
%create new variables that start from t = 0 and only contain event data
newtime=time(starttime:length(time))-time(starttime);

newvoltage=voltage(starttime:length(time));

Tstart = time(starttime);

Vstart = voltage(starttime);


% output = [tout,vout];
output = [newtime,newvoltage];
end

