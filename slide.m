function [output] = slide(time,voltage)
%outputs the maximum slope and its position in the matrix
%time is x, voltage is y
load lab2part1variables.mat

for i = 1:1:(length(time)-51)
    limit = (i+51);
    tmask = time(i:limit);
    vmask = voltage(i:limit);
    fit = polyfit(tmask,vmask,1);
    m = fit(1);
    slopes(i) = m; %an array of slopes at each point in the line
    posslopes = abs(slopes);
    [~,pos] = max(posslopes);
    maxslope = slopes(pos);
    
    newtime = time(pos:length(time));
    time = time-time(pos);
    newvoltage = voltage(pos:length(voltage));
    
    tcv = (voltage*1000-betaHat(1))/betaHat(2); %°C newvoltage for trimmed, voltage for untrimmed
    
    output = [time,tcv];
end



end

