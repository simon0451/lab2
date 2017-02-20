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
%disp(betaHat);
% plot the best fit line
xx = linspace(0,100);
yy = betaHat(1) + betaHat(2)*xx; %betaHat(1) is the Y-intercept, and betaHat(2) is the slope
% plot the points (data) for which we found the best fit
m = num2str(betaHat(2),3);
b = num2str(betaHat(1),3);
txt = strcat('y=',m,'x+',b);
txt1 = strcat('Units: (mV)=(mV/°C)*(°C)+(mV)');

figure(1)
plot(ThermistorTemperature,ThermocoupleVoltage,'o',xx,yy,'--')
title('Thermocouple Measurement Response')
ylabel('Thermocouple Output Voltage (mV)')
xlabel('Temperature From Thermistor (°C)')
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
ThermocoupleTemperature = (ThermocoupleVoltage-betaHat(1))/betaHat(2); %°C
Part3BF = polyfit(ThermistorTemperature,ThermocoupleTemperature,1);
p3bfyvalues = Part3BF(1)*ThermistorTemperature+Part3BF(2);

Yc = p3bfyvalues; %The value of y predicted by the polynomial equation for a given value of x
tvp = 2.262; %For N = 10, 95% confidence
yiyci = (ThermocoupleTemperature-p3bfyvalues).^2;
sumyiyci = sum(yiyci);
Syx = (sumyiyci/(length(ThermocoupleTemperature)-1))^.5; %standard error of the fit
SampleMeanValue = (sum(ThermistorTemperature))/length(ThermistorTemperature);

for i = 1:1:length(ThermistorTemperature)
    unsummedDen(i) = (ThermistorTemperature(i)-SampleMeanValue)^2;
end
Den = sum(unsummedDen);

CIofFitPOS = Yc+tvp.*Syx.*(1./length(ThermocoupleTemperature)+((ThermistorTemperature-SampleMeanValue).^2./(Den))).^.5;
CIofFitNEG = Yc-tvp.*Syx.*(1./length(ThermocoupleTemperature)+((ThermistorTemperature-SampleMeanValue).^2./(Den))).^.5;

CIofMeasurementPOS = Yc+tvp.*Syx.*(1+1./length(ThermocoupleTemperature)+((ThermistorTemperature-SampleMeanValue).^2./(Den))).^.5;
CIofMeasurementNEG = Yc-tvp.*Syx.*(1+1./length(ThermocoupleTemperature)+((ThermistorTemperature-SampleMeanValue).^2./(Den))).^.5;


figure(2)
plot(ThermistorTemperature,ThermocoupleTemperature,'o',ThermistorTemperature,p3bfyvalues,'--',ThermistorTemperature,CIofFitPOS,'y--',ThermistorTemperature,CIofMeasurementNEG,'r--',ThermistorTemperature,CIofFitNEG,'y--',ThermistorTemperature,CIofMeasurementPOS,'r--')
title('Thermocouple Response vs. Thermistor Response')
ylabel('Temperature From Thermocouple (°C)')
xlabel('Temperature From Thermistor (°C)')
grid on
xmin = -5;
xmax = 105;
ymin = -5;
ymax = 106;
axis ([xmin xmax ymin ymax])
legend('Thermocouple Output','Linear Least Squares Fit','Confidence Interval of Fit','Confidence Interval of Measurement','location','southeast')

figure(3)
plot(ThermistorTemperature,ThermocoupleTemperature,'o',ThermistorTemperature,p3bfyvalues,'--',ThermistorTemperature,CIofFitPOS,'y--',ThermistorTemperature,CIofMeasurementNEG,'r--',ThermistorTemperature,CIofFitNEG,'y--',ThermistorTemperature,CIofMeasurementPOS,'r--')
title('Zoomed-In Thermocouple Response vs. Thermistor Response')
ylabel('Temperature From Thermocouple (°C)')
xlabel('Temperature From Thermistor (°C)')
grid on
xmin = 30;
xmax = 50;
ymin = 30;
ymax = 50;
axis ([xmin xmax ymin ymax])
legend('Thermocouple Output','Linear Least Squares Fit','Confidence Interval of Fit','Confidence Interval of Measurement','location','southeast')

%% Static Calibration Part 4
Thermo25Temperature = (TCV-betaHat(1))/betaHat(2); %°C
Tbar = (sum(Thermo25Temperature))/length(Thermo25Temperature); %sample mean value
StandardDeviation25 = std(Thermo25Temperature);
N = length(Thermo25Temperature);
v = N-1;
tvp25 = 2.067; %95% confidence
AM = sum(Thermo25Temperature)/length(Thermo25Temperature); %arithmetic mean
for i = 1:1:length(Thermo25Temperature)
    sxcomp(i) = (Thermo25Temperature(i)-AM)^2;
end
compsum = sum(sxcomp);
Sx = ((1/v)*compsum)^.5;
Sxbar = Sx/((N)^.5);
XiSPOS = AM+tvp25*Sx; %positive 95% confidence limit of measurement (outer lines)
XiSNEG = AM-tvp25*Sx;
XiPOS = AM+tvp25*Sxbar; %positive 95% confidence limit (true mean value - inner lines)
XiNEG = AM-tvp25*Sxbar;

%% Static Calibration Part 5

for i = 1:1:25
    zeroC(i) = 0;
end
zeroC = zeroC';

figure(4)
plot(ThermistorTemperature,ThermocoupleTemperature,'o',ThermistorTemperature,p3bfyvalues,'--',ThermistorTemperature,CIofFitPOS,'y--',ThermistorTemperature,CIofMeasurementNEG,'r--',zeroC,Thermo25Temperature,'bo',ThermistorTemperature,CIofFitNEG,'y--',ThermistorTemperature,CIofMeasurementPOS,'r--')
title('Thermocouple Response vs. Thermistor Response')
ylabel('Temperature From Thermocouple (°C)')
xlabel('Temperature From Thermistor (°C)')
grid on
xmin = -2;
xmax = 2;
ymin = -2;
ymax = 2;
axis ([xmin xmax ymin ymax])
legend('Thermocouple Output','Linear Least Squares Fit','Confidence Interval of Fit','Confidence Interval of Measurement','location','southeast')


