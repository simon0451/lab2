function [output] = pros(time,voltage,startingtime)
offset = 0;

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

