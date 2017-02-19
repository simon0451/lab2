%% Header

%Simon Popecki
%19 February 2017
%ME 646
%Lab 2

%REQUIRED FILES:
%lab2.mat

%% Static Calibration Part 2
clear all, close all;
load lab2.mat
%TRC has the following units: deg. C, KOhms, mV
%Column A is temperature of the bath, column B is the resistance of the
%thermistor, and column C is the voltage of the thermocouple.
Ro = 9.64788; %kOhms
B = 3617.58; %units????
To = 298.15; %K
ThermistorResistance = TRC(:,2); %kOhms
ThermistorInverseT = (1/To)+(1/B).*log(ThermistorResistance./Ro);
ThermistorTemperatureKelvin = 1./ThermistorInverseT; %K
ThermistorTemperature = ThermistorTemperatureKelvin-273.15; %C
ThermocoupleVoltage = TRC(:,3); %mV

%Least squares fit
%Example code source: https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)
input = [ThermistorTemperature,ThermocoupleVoltage]; %input line
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
txt = strcat('y=',m,'x+',b);
txt1 = strcat('Units: (mV)=(mV/�C)*(�C)+(mV)');

figure(1)
plot(ThermistorTemperature,ThermocoupleVoltage,'o',xx,yy,'--')
title('Thermocouple Measurement Response')
ylabel('Thermocouple Output Voltage (mV)')
xlabel('Temperature From Thermistor (�C)')
grid on
xmin = -5;
xmax = 105;
ymin = -100;
ymax = 1100;
axis ([xmin xmax ymin ymax])
text(.55*xmax,.18*ymax,txt)
text(.55*xmax,.13*ymax,txt1)
legend('Thermocouple Output','Linear Least Squares Fit','location','southeast')

%% Static Calibration Part 3

%% Static Calibration Part 4








