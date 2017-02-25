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

bareboilicetime = xlsread('Michalak_Popecki_Rose.xlsx',8,'a:a');
bareboilicevoltage = xlsread('Michalak_Popecki_Rose.xlsx',8,'B9:B5008');

bareiceboiltime = xlsread('Michalak_Popecki_Rose.xlsx',9,'a:a');
bareiceboilvoltage = xlsread('Michalak_Popecki_Rose.xlsx',9,'B9:B12008');

%REFERENCE PLOTS:
% figure(5)
% plot(steelboilicetime,steelboilicevoltage,alumboilicetime,alumboilicevoltage,bareboilicetime,bareboilicevoltage)
% title('(Un-processed Data) Thermocouples - Boiling Water to Ice Water')
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple')
% grid on
% 
% figure(6)
% plot(steeliceboiltime,steeliceboilvoltage,alumiceboiltime,alumiceboilvoltage,bareiceboiltime,bareiceboilvoltage)
% title('(Un-processed Data) Thermocouples - Ice Water to Boiling Water')
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple','location','southeast')
% grid on

%determining start of data position using the 5-sigma method, smoothing

%smoothing data:
%the 51 represents the mask width
steelboilicevoltage = smooth(steelboilicevoltage,51);
alumboilicevoltage = smooth(alumboilicevoltage,51);
bareboilicevoltage = smooth(bareboilicevoltage,51);

steeliceboilvoltage = smooth(steeliceboilvoltage,51);
alumiceboilvoltage = smooth(alumiceboilvoltage,51);
bareiceboilvoltage = smooth(bareiceboilvoltage,51);

%boiling water to ice water - METHOD 1
%Using the tuning factor: using the wrong tuning factor will either throw
%an error response or result in the data not being started at the proper
%time (usually the idle time in the beginning is not cut off like it should
%be). The tuning factor should be adjusted to the poin where the input
%function when drawn on a plot, "snaps" to the starting point.
steelboilicearray = pros(steelboilicetime,steelboilicevoltage,1); %outputs [time,temperature, start time tuning factor] of the input using method 1
alumboilicearray = pros(alumboilicetime,alumboilicevoltage,.5);
bareboilicearray = pros(bareboilicetime,bareboilicevoltage,0);
%ice water to boiling water
steeliceboilarray = pros(steeliceboiltime,steeliceboilvoltage,.5);
alumiceboilarray = pros(alumiceboiltime,alumiceboilvoltage,.5);
bareiceboilarray = pros(bareiceboiltime,bareiceboilvoltage,1.41);

%boiling water to icewater - METHOD 2
steelboilicearray2 = slide(steelboilicetime,steelboilicevoltage);
alumboilicearray2 = slide(alumboilicetime,alumboilicevoltage);
bareboilicearray2 = slide(bareboilicetime,bareboilicevoltage);

%ice water to boiling water - METHOD 2
steeliceboilarray2 = slide(steeliceboiltime,steeliceboilvoltage);
alumiceboilarray2 = slide(alumiceboiltime,alumiceboilvoltage);
bareiceboilarray2 = slide(bareiceboiltime,bareiceboilvoltage);

%finding Tfinal for the bare wire thermocouples going from ice water to
%boiling water
%Time measurements are in .001 second intervals
%averaging the last 2 seconds = last 2,000 measurements of the array - 1.5
%s = 1,500 measurements
twosec = 1500;
bareiceboilT = bareiceboilarray(:,2);
lbareiceboilT = length(bareiceboilT); %some number like 3612 - the length of the vector
bareiceboildatastart = lbareiceboilT-twosec; %the position in the array where we begin looking at data
bareiceboilrange = bareiceboilT(bareiceboildatastart:lbareiceboilT);
bareiceboilTfinal = mean(bareiceboilrange); %the average temperature of the boiling water bath, celcius

%for boiling water to ice water
bareboiliceT =bareboilicearray(:,2);
lbareboiliceT = length(bareboiliceT);
bareboilicedatastart = lbareboiliceT-twosec;
bareboilicerange = bareboiliceT(bareboilicedatastart:lbareboiliceT);
bareboiliceTfinal = mean(bareboilicerange); %degrees celcius

%embedded thermocouoples final temperatures
%using the final value instead of averaging

%for ice water to boiling water
steeliceboilT = steeliceboilarray(:,2);
steeliceboilTfinal = steeliceboilT(end);

alumiceboilT = alumiceboilarray(:,2);
alumiceboilTfinal = alumiceboilT(end);

%for boiling water to ice water
steelboiliceT = steelboilicearray(:,2);
steelboiliceTfinal = steelboiliceT(end);

alumboiliceT = alumboilicearray(:,2);
alumboiliceTfinal = alumboiliceT(end);

%plottting boiling water to ice water transfer
%NOTE! Individual points are NOT being plotted - they form a thick line and
%it looks terrible.
figure(7)
plot(steelboilicearray(:,1),steelboilicearray(:,2),alumboilicearray(:,1),alumboilicearray(:,2),bareboilicearray(:,1),bareboilicearray(:,2))
title('(Method 1) Thermocouples - Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (C)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple')
grid on

%plotting ice water to boiling water transfer
figure(8)
plot(steeliceboilarray(:,1),steeliceboilarray(:,2),alumiceboilarray(:,1),alumiceboilarray(:,2),bareiceboilarray(:,1),bareiceboilarray(:,2))
title('(Method 1) Thermocouples - Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (C)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple','location','southeast')
grid on

figure(9)
plot(steelboilicearray2(:,1),steelboilicearray2(:,2),alumboilicearray2(:,1),alumboilicearray2(:,2),bareboilicearray2(:,1),bareboilicearray2(:,2))
title('(Method 2) Thermocouples - Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (C)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple')
grid on

%plotting ice water to boiling water transfer
figure(10)
plot(steeliceboilarray2(:,1),steeliceboilarray2(:,2),alumiceboilarray2(:,1),alumiceboilarray2(:,2),bareiceboilarray2(:,1),bareiceboilarray2(:,2))
title('(Method 2) Thermocouples - Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (C)')
legend('Steel Embedded Thermocouple','Aluminum Embedded Thermocouple','Bare Wire Thermocouple','location','southeast')
grid on

%% Dynamic Calibration Part 2
load d1variables.mat
%solving for the time constants

%ice water to boiling water
%the timecon function takes an input in C but works in K, the answer is
%normalized, so the temperature units do not matter - it is to avoid a
%complex number answer

%The resutls are the log natural of Gamma
lnGammasteeliceboil = timecon(steeliceboilarray(:,2),steeliceboilTfinal); %(the temperature function in C, the final temperature)
lnGammasteeliceboil2 = timecon(steeliceboilarray2(:,2),steeliceboilTfinal);

lnGammaalumiceboil = timecon(alumiceboilarray(:,2),alumiceboilTfinal);
lnGammaalumiceboil2 = timecon(alumiceboilarray2(:,2),alumiceboilTfinal);

lnGammabareiceboil = timecon(bareiceboilarray(:,2),bareiceboilTfinal);
lnGammabareiceboil2 = timecon(bareiceboilarray2(:,2),bareiceboilTfinal);

%the same results for the transition from boiling water to ice water
lnGammasteelboilice = timecon(steelboilicearray(:,2),steelboiliceTfinal); %%%%%%%%%%%%%%%%%%%%%%
lnGammasteelboilice2 = timecon(steelboilicearray2(:,2),steelboiliceTfinal);

lnGammaalumboilice = timecon(alumboilicearray(:,2),alumboiliceTfinal);
lnGammaalumboilice2 = timecon(alumboilicearray2(:,2),alumboiliceTfinal);

lnGammabareboilice = timecon(bareboilicearray(:,2),bareboiliceTfinal);
lnGammabareboilice2 = timecon(bareboilicearray2(:,2),bareboiliceTfinal);

tausteelboilice = -(steelboilicearray(:,1))./(lnGammasteelboilice);

figure(21)
subplot(3,1,1)
plot(steelboilicearray(:,1),(lnGammasteelboilice))
title('Method 1, Thermocouples - Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Natural log of Gamma','location','northeast')
xmin = -5;
xmax = 50;
ymin = -5;
ymax = 1;
axis ([xmin xmax ymin ymax])
grid on




%% Dynamic Calibration Part 5
load lab2part1variables.mat

bareiceairtime = xlsread('Michalak_Popecki_Rose.xlsx',5,'a:a'); %this is for bareiceair3 - three samples were taken and this one has the best data
bareiceairvoltage = xlsread('Michalak_Popecki_Rose.xlsx',5,'B9:B12008');
%a very noisy signal...

bareiceairvoltage = smooth(bareiceairvoltage,1001);
bareiceairarray = pros(bareiceairtime,bareiceairvoltage,2.19);
%bareiceairtemperature = ((bareiceairarray(:,2))-betaHat(1))/betaHat(2); %betahat 2 is the slope

rt = 26; %°C, from my lab notebook
for i = 1:1:length(bareiceairarray(:,1))
    roomtemp(i) = rt;
end
    
figure(51)
plot(bareiceairarray(:,1),bareiceairarray(:,2),bareiceairarray(:,1),roomtemp,'--')
title('Bare Wire Thermocouple - Ice Water to Air')
xlabel('Time (s)')
ylabel('Temperature (C)')
legend('Bare Wire Thermocouple Temperature','Room Temperature','location','southeast')
grid on



























