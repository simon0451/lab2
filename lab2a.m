%% Header

%Simon Popecki
%19 February 2017
%ME 646
%Lab 2

%% Static Calibration Part 2
clear all, close all;
load lab2.mat
%TRC has the following units: deg. C, KOhms, mV
%Column A is temperature of the bath, column B is the resistance of the
%thermistor, and column C is the voltage of the thermocouple.

BathTemperature = TRC(:,1); %C
ThermocoupleVoltage = TRC(:,3); %mV

%Least squares fit
%Example code source: https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)

input = [BathTemperature,ThermocoupleVoltage]; %input line
pts = length(input);             % number of points
X = [ones(pts,1), input(:,1)];   % forming X of X beta = y
y = input(:,2);                % forming y of X beta = y
betaHat = (X' * X) \ X' * y;   % computing projection of matrix X on y, giving beta
% display best fit parameters
disp(betaHat);
% plot the best fit line
xx = linspace(0,100);
yy = betaHat(1) + betaHat(2)*xx; %betaHat(1) is the Y-intercept, and betaHat(2) is the slope
% plot the points (data) for which we found the best fit
m = num2str(betaHat(2),3);
b = num2str(betaHat(1),3);
txt = strcat('y (mV) = ',m,' (mV/�C) x (�C) +',b,' (mV)');

figure(1)
plot(BathTemperature,ThermocoupleVoltage,'o',xx,yy,'--')
title('Thermocouple Measurement Response')
ylabel('Voltage (mV)')
xlabel('Temperature (�C)')
grid on
xmin = -5;
xmax = 105;
ymin = -100;
ymax = 1100;
axis ([xmin xmax ymin ymax])
text(.45*xmax,.13*ymax,txt)
legend('Thermocouple Output','Linear Least Squares Fit','location','southeast')





