function [slopes] = slide(time,voltage)
%time is x, voltage is y


for i = 1:1:(length(time)-51)
    limit = (i+51);
    tmask = time(i:limit);
    vmask = voltage(i:limit);
    fit = polyfit(tmask,vmask,1);
    m = fit(1);
    slopes(i) = m; %an array of slopes at each point in the line
end



end

