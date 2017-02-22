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
tvp25 = 2.067; %95% confidence, from table
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

%% Dynamic Calibration Part 1
%Determining the time that the thermocouples transition to the new bath
%cleaning up data using two methods

clear all;
%time is in seconds
%voltage is in volts

%loading in the data from the excel file
%time is in seconds
%voltage is in volts
steelboilicetime = xlsread('Michalak_Popecki_Rose.xlsx',1,'a:a'); %time is the first column on each measurement, can be different on different measurements
steelboilicevoltage = xlsread('Michalak_Popecki_Rose.xlsx',1,'B9:B5008');

alumboilicetime = xlsread('Michalak_Popecki_Rose.xlsx',2,'a:a');
alumboilicevoltage = xlsread('Michalak_Popecki_Rose.xlsx',2,'B9:B5008');

steeliceboiltime = xlsread('Michalak_Popecki_Rose.xlsx',3,'a:a');
steeliceboilvoltage = xlsread('Michalak_Popecki_Rose.xlsx',3,'B9:B5008');

alumiceboiltime = xlsread('Michalak_Popecki_Rose.xlsx',4,'a:a');
alumiceboilvoltage = xlsread('Michalak_Popecki_Rose.xlsx',4,'B9:B5008');

bareiceairtime = xlsread('Michalak_Popecki_Rose.xlsx',5,'a:a'); %this is for bareiceair3 - three samples were taken and this one has the best data
bareiceairvoltage = xlsread('Michalak_Popecki_Rose.xlsx',5,'B9:B12008');

bareboilicetime = xlsread('Michalak_Popecki_Rose.xlsx',8,'a:a');
bareboilicevoltage = xlsread('Michalak_Popecki_Rose.xlsx',8,'B9:B5008');

bareiceboiltime = xlsread('Michalak_Popecki_Rose.xlsx',9,'a:a');
bareiceboilvoltage = xlsread('Michalak_Popecki_Rose.xlsx',9,'B9:B12008');

figure(5)
plot(steelboilicetime,steelboilicevoltage,alumboilicetime,alumboilicevoltage,bareboilicetime,bareboilicevoltage)
title('(Un-processed Data) Thermocouples - Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple')
grid on

figure(6)
plot(steeliceboiltime,steeliceboilvoltage,alumiceboiltime,alumiceboilvoltage,bareiceboiltime,bareiceboilvoltage)
title('(Un-processed Data) Thermocouples - Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple','location','southeast')
grid on

figure(7)
plot(bareiceairtime,bareiceairvoltage)
title('(Un-processed Data) Bare Wire Thermocouple - Ice Water to Air')
xlabel('Time (s)')
ylabel('Voltage (V)')
grid on

%determining start of data position using the 5-sigma method, smoothing

%boiling water to ice water
%Using the tuning factor: using the wrong tuning factor will either throw
%an error response or result in the data not being started at the proper
%time (usually the idle time in the beginning is not cut off like it should
%be). The tuning factor should be adjusted to the poin where the input
%function when drawn on a plot, "snaps" to the starting point.
steelboilicearray = pros(steelboilicetime,steelboilicevoltage,1); %outputs [time,voltage, start time tuning factor] of the input using method 1
alumboilicearray = pros(alumboilicetime,alumboilicevoltage,.5);
bareboilicearray = pros(bareboilicetime,bareboilicevoltage,0);

%ice water to boiling water
steeliceboilarray = pros(steeliceboiltime,steeliceboilvoltage,.5);
alumiceboilarray = pros(alumiceboiltime,alumiceboilvoltage,.5);
bareiceboilarray = pros(bareiceboiltime,bareiceboilvoltage,1.41);

%plottting boiling water to ice water transfer

%NOTE! Individual points are NOT being plotted - they form a thick line and
%it looks terrible.
figure(8)
plot(steelboilicearray(:,1),steelboilicearray(:,2),alumboilicearray(:,1),alumboilicearray(:,2),bareboilicearray(:,1),bareboilicearray(:,2))
title('(Method 1) Thermocouples - Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple')
grid on

%plotting ice water to boiling water transfer
figure(9)
plot(steeliceboilarray(:,1),steeliceboilarray(:,2),alumiceboilarray(:,1),alumiceboilarray(:,2),bareiceboilarray(:,1),bareiceboilarray(:,2))
title('(Method 1) Thermocouples - Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Voltage (V)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple','location','southeast')
grid on































