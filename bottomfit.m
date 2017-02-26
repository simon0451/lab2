function [output] = bottomfit(middlefit)
%This provides an array to plot the residuals vs. time

%the residuals are the difference between the data and the prediction
%this function will use the results of "middlefit.m" as an input

datay = middlefit(:,2); %the second column of the middlefit output array
predy = middlefit(:,4); %the predicted y values
xdata = middlefit(:,1); %remember that the xvalues are the same for both arrays

residuals = datay - predy; %an array of differences between the actual data and the prediction

output = [xdata,residuals];



end

