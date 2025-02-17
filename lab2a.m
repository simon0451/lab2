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
txt = strcat('V (mV) =',m,' (mV/�C) T (�C)+',b,' (mV)');

figure(1)
plot(ThermistorTemperature,ThermocoupleVoltage,'o',xx,yy,'--')
ylabel('Thermocouple Output Voltage (mV)')
xlabel('Temperature From Thermistor (�C)')
grid on
xmin = -5;
xmax = 105;
ymin = -100;
ymax = 1100;
axis ([xmin xmax ymin ymax])
text(.45*xmax,.18*ymax,txt)
legend('Thermocouple Output','Linear Least Squares Fit','location','southeast')

%% Static Calibration Part 3
ThermocoupleTemperature = (ThermocoupleVoltage-betaHat(1))/betaHat(2); %�C
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
ylabel('Temperature From Thermocouple (�C)')
xlabel('Temperature From Thermistor (�C)')
grid on
xmin = -5;
xmax = 105;
ymin = -5;
ymax = 106;
axis ([xmin xmax ymin ymax])
legend('Thermocouple Output','Linear Least Squares Fit','Confidence Interval of Fit','Confidence Interval of Measurement','location','southeast')

% figure(3)
% plot(ThermistorTemperature,ThermocoupleTemperature,'o',ThermistorTemperature,p3bfyvalues,'--',ThermistorTemperature,CIofFitPOS,'y--',ThermistorTemperature,CIofMeasurementNEG,'r--',ThermistorTemperature,CIofFitNEG,'y--',ThermistorTemperature,CIofMeasurementPOS,'r--')
% title('Zoomed-In Thermocouple Response vs. Thermistor Response')
% ylabel('Temperature From Thermocouple (�C)')
% xlabel('Temperature From Thermistor (�C)')
% grid on
% xmin = 30;
% xmax = 50;
% ymin = 30;
% ymax = 50;
% axis ([xmin xmax ymin ymax])
% legend('Thermocouple Output','Linear Least Squares Fit','Confidence Interval of Fit','Confidence Interval of Measurement','location','southeast')

%% Static Calibration Part 4
Thermo25Temperature = (TCV-betaHat(1))/betaHat(2); %�C
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
plot(ThermistorTemperature,ThermocoupleTemperature,'o',ThermistorTemperature,p3bfyvalues,'--',ThermistorTemperature,CIofFitPOS,'b--',ThermistorTemperature,CIofMeasurementNEG,'r--',zeroC,Thermo25Temperature,'bo',ThermistorTemperature,CIofFitNEG,'b--',ThermistorTemperature,CIofMeasurementPOS,'r--')
title('Thermocouple Response vs. Thermistor Response')
ylabel('Temperature From Thermocouple (�C)')
xlabel('Temperature From Thermistor (�C)')
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

%% STEEL FROM BOILING WATER TO ICE WATER
partisbi = gammafit(steelboilicearray(:,1),steelboilicearray(:,2));
partiisbi = middlefit(steelboilicearray(:,1),steelboilicearray(:,2));
partiiisbi = bottomfit(partiisbi);
D3sbi = p2(steelboilicearray(:,1),steelboilicearray(:,2));

figure(31)
subplot(2,1,1)
plot(D3sbi(:,1),D3sbi(:,2),D3sbi(:,3),D3sbi(:,4))
title('5 \sigma Method, Steel Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction')
grid on
axis([-3 40 -10 110])
subplot(2,1,2)
plot(D3sbi(:,5),D3sbi(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 40 -20 20])

figure(21)
subplot(3,1,1)
plot(partisbi(:,1),partisbi(:,2),partisbi(:,3),partisbi(:,4))
title('5 \sigma Method, Steel Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 30 -inf 1])
subplot(3,1,2)
plot(partiisbi(:,1),partiisbi(:,2),partiisbi(:,3),partiisbi(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 30 -10 110])
grid on
subplot(3,1,3)
plot(partiiisbi(:,1),partiiisbi(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 30 -20 20])

%USING METHOD 2
partisbi2 = gammafit(steelboilicearray2(:,1),steelboilicearray2(:,2));
partiisbi2 = middlefit(steelboilicearray2(:,1),steelboilicearray2(:,2));
partiiisbi2 = bottomfit(partiisbi2);
D3sbi2 = p2(steelboilicearray2(:,1),steelboilicearray2(:,2));

figure(32)
subplot(2,1,1)
plot(D3sbi2(:,1),D3sbi2(:,2),D3sbi2(:,3),D3sbi2(:,4))
title('Max Slope Method, Steel Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction')
grid on
axis([-3 40 -10 110])
subplot(2,1,2)
plot(D3sbi2(:,5),D3sbi2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 40 -20 20])

figure(22)
subplot(3,1,1)
plot(partisbi2(:,1),partisbi2(:,2),partisbi2(:,3),partisbi2(:,4))
title('Max Slope Method, Steel Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 30 -inf 1])
subplot(3,1,2)
plot(partiisbi2(:,1),partiisbi2(:,2),partiisbi2(:,3),partiisbi2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 30 -10 110])
grid on
subplot(3,1,3)
plot(partiiisbi2(:,1),partiiisbi2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 30 -20 20])

%% STEEL FROM ICE WATER TO BOILING WATER
partisib = gammafit(steeliceboilarray(:,1),steeliceboilarray(:,2));
partiisib = middlefit(steeliceboilarray(:,1),steeliceboilarray(:,2));
partiiisib = bottomfit(partiisib);
D3sib = p2(steeliceboilarray(:,1),steeliceboilarray(:,2));

figure(33)
subplot(2,1,1)
plot(D3sib(:,1),D3sib(:,2),D3sib(:,3),D3sib(:,4))
title('5 \sigma Method, Steel Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3sib(:,5),D3sib(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(23)
subplot(3,1,1)
plot(partisib(:,1),partisib(:,2),partisib(:,3),partisib(:,4))
title('5 \sigma Method, Steel Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiisib(:,1),partiisib(:,2),partiisib(:,3),partiisib(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiisib(:,1),partiiisib(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])

partisib2 = gammafit(steeliceboilarray2(:,1),steeliceboilarray2(:,2));
partiisib2 = middlefit(steeliceboilarray2(:,1),steeliceboilarray2(:,2));
partiiisib2 = bottomfit(partiisib2);
D3sib2 = p2(steeliceboilarray2(:,1),steeliceboilarray2(:,2));

figure(34)
subplot(2,1,1)
plot(D3sib2(:,1),D3sib2(:,2),D3sib2(:,3),D3sib2(:,4))
title('Max Slope Method, Steel Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3sib2(:,5),D3sib2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(24)
subplot(3,1,1)
plot(partisib2(:,1),partisib2(:,2),partisib2(:,3),partisib2(:,4))
title('Max Slope Method, Steel Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiisib2(:,1),partiisib2(:,2),partiisib2(:,3),partiisib2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiisib2(:,1),partiiisib2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])

%% ALUMINUM BOILING WATER TO ICE WATER

partiabi = gammafit(alumboilicearray(:,1),alumboilicearray(:,2));
partiiabi = middlefit(alumboilicearray(:,1),alumboilicearray(:,2));
partiiiabi = bottomfit(partiiabi);
D3abi = p2(alumboilicearray(:,1),alumboilicearray(:,2));

figure(35)
subplot(2,1,1)
plot(D3abi(:,1),D3abi(:,2),D3abi(:,3),D3abi(:,4))
title('5 \sigma Method, Aluminum Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3abi(:,5),D3abi(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(25)
subplot(3,1,1)
plot(partiabi(:,1),partiabi(:,2),partiabi(:,3),partiabi(:,4))
title('5 \sigma Method, Aluminum Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiiabi(:,1),partiiabi(:,2),partiiabi(:,3),partiiabi(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiiabi(:,1),partiiiabi(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])

partiabi2 = gammafit(alumboilicearray2(:,1),alumboilicearray2(:,2));
partiiabi2 = middlefit(alumboilicearray2(:,1),alumboilicearray2(:,2));
partiiiabi2 = bottomfit(partiiabi2);
D3abi2 = p2(alumboilicearray2(:,1),alumboilicearray2(:,2));

figure(36)
subplot(2,1,1)
plot(D3abi2(:,1),D3abi2(:,2),D3abi2(:,3),D3abi2(:,4))
title('Max Slope Method, Aluminum Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3abi2(:,5),D3abi2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(26)
subplot(3,1,1)
plot(partiabi2(:,1),partiabi2(:,2),partiabi2(:,3),partiabi2(:,4))
title('Max Slope Method, Aluminum Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiiabi2(:,1),partiiabi2(:,2),partiiabi2(:,3),partiiabi2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiiabi2(:,1),partiiiabi2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])


%% ALUMINUM ICE WATER TO BOILING WATER
partiaib = gammafit(alumiceboilarray(:,1),alumiceboilarray(:,2));
partiiaib = middlefit(alumiceboilarray(:,1),alumiceboilarray(:,2));
partiiiaib = bottomfit(partiiaib);
D3aib = p2(alumiceboilarray(:,1),alumiceboilarray(:,2));

figure(37)
subplot(2,1,1)
plot(D3aib(:,1),D3aib(:,2),D3aib(:,3),D3aib(:,4))
title('5 \sigma Method, Aluminum Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3aib(:,5),D3aib(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(27)
subplot(3,1,1)
plot(partiaib(:,1),partiaib(:,2),partiaib(:,3),partiaib(:,4))
title('5 \sigma Method, Aluminum Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiiaib(:,1),partiiaib(:,2),partiiaib(:,3),partiiaib(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiiaib(:,1),partiiiaib(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])

partiaib2 = gammafit(alumiceboilarray2(:,1),alumiceboilarray2(:,2));
partiiaib2 = middlefit(alumiceboilarray2(:,1),alumiceboilarray2(:,2));
partiiiaib2 = bottomfit(partiiaib2);
D3aib2 = p2(alumiceboilarray2(:,1),alumiceboilarray2(:,2));

figure(38)
subplot(2,1,1)
plot(D3aib2(:,1),D3aib2(:,2),D3aib2(:,3),D3aib2(:,4))
title('Max Slope Method, Aluminum Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-3 20 -10 110])
subplot(2,1,2)
plot(D3aib2(:,5),D3aib2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-3 20 -20 20])

figure(28)
subplot(3,1,1)
plot(partiaib2(:,1),partiaib2(:,2),partiaib2(:,3),partiaib2(:,4))
title('Max Slope Method, Aluminum Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])
subplot(3,1,2)
plot(partiiaib2(:,1),partiiaib2(:,2),partiiaib2(:,3),partiiaib2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-5 20 -10 110])
grid on
subplot(3,1,3)
plot(partiiiaib2(:,1),partiiiaib2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-5 20 -20 20])

%% BARE WIRE BOILING WATER TO ICE WATER


partibbi = gammafit(bareboilicearray(:,1),bareboilicearray(:,2));
partiibbi = middlefit(bareboilicearray(:,1),bareboilicearray(:,2));
partiiibbi = bottomfit(partiibbi);
D3bbi = p2(bareboilicearray(:,1),bareboilicearray(:,2));

figure(39)
subplot(2,1,1)
plot(D3bbi(:,1),D3bbi(:,2),D3bbi(:,3),D3bbi(:,4))
title('5 \sigma Method, Bare Wire Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction')
grid on
axis([-1 1 -10 110])
subplot(2,1,2)
plot(D3bbi(:,5),D3bbi(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-1 1 -20 20])

figure(29)
subplot(3,1,1)
plot(partibbi(:,1),partibbi(:,2),partibbi(:,3),partibbi(:,4))
title('5 \sigma Method, Bare Wire Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-1 1 -inf 1])
subplot(3,1,2)
plot(partiibbi(:,1),partiibbi(:,2),partiibbi(:,3),partiibbi(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-1 1 -10 110])
grid on
subplot(3,1,3)
plot(partiiibbi(:,1),partiiibbi(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-1 1 -20 20])

%USING METHOD 2
partibbi2 = gammafit(bareboilicearray2(:,1),bareboilicearray2(:,2));
partiibbi2 = middlefit(bareboilicearray2(:,1),bareboilicearray2(:,2));
partiiibbi2 = bottomfit(partiibbi2);
D3bbi2 = p2(bareboilicearray2(:,1),bareboilicearray2(:,2));

figure(310)
subplot(2,1,1)
plot(D3bbi2(:,1),D3bbi2(:,2),D3bbi2(:,3),D3bbi2(:,4))
title('Max Slope Method, Bare Wire Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction')
grid on
axis([-1 1 -10 110])
subplot(2,1,2)
plot(D3bbi2(:,5),D3bbi2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-1 1 -20 20])

figure(210)
subplot(3,1,1)
plot(partibbi2(:,1),partibbi2(:,2),partibbi2(:,3),partibbi2(:,4))
title('Max Slope Method, Bare Wire Embedded Thermocouple, Boiling Water to Ice Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-1 1 -inf 1])
subplot(3,1,2)
plot(partiibbi2(:,1),partiibbi2(:,2),partiibbi2(:,3),partiibbi2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-1 1 -10 110])
grid on
subplot(3,1,3)
plot(partiiibbi2(:,1),partiiibbi2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-1 1 -20 20])

%% BARE WIRE ICE WATER TO BOILING WATER

partibib = gammafit(bareiceboilarray(:,1),bareiceboilarray(:,2));
partiibib = middlefit(bareiceboilarray(:,1),bareiceboilarray(:,2));
partiiibib = bottomfit(partiibib);
D3bib = p2(bareiceboilarray(:,1),bareiceboilarray(:,2));

figure(311)
subplot(2,1,1)
plot(D3bib(:,1),D3bib(:,2),D3bib(:,3),D3bib(:,4))
title('5 \sigma Method, Bare Wire Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-1 1 -10 110])
subplot(2,1,2)
plot(D3bib(:,5),D3bib(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-1 1 -20 20])

figure(211)
subplot(3,1,1)
plot(partibib(:,1),partibib(:,2),partibib(:,3),partibib(:,4))
title('5 \sigma Method, Bare Wire Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-1 1 -inf 1])
subplot(3,1,2)
plot(partiibib(:,1),partiibib(:,2),partiibib(:,3),partiibib(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-1 1 -10 110])
grid on
subplot(3,1,3)
plot(partiiibib(:,1),partiiibib(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-1 1 -20 20])

partibib2 = gammafit(bareiceboilarray2(:,1),bareiceboilarray2(:,2));
partiibib2 = middlefit(bareiceboilarray2(:,1),bareiceboilarray2(:,2));
partiiibib2 = bottomfit(partiibib2);
D3bib2 = p2(bareiceboilarray2(:,1),bareiceboilarray2(:,2));

figure(312)
subplot(2,1,1)
plot(D3bib2(:,1),D3bib2(:,2),D3bib2(:,3),D3bib2(:,4))
title('Max Slope Method, Bare Wire Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','southeast')
grid on
axis([-1 1 -10 110])
subplot(2,1,2)
plot(D3bib2(:,5),D3bib2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
legend('Data Minus Prediction','location','southeast')
grid on
axis([-1 1 -20 20])

figure(212)
subplot(3,1,1)
plot(partibib2(:,1),partibib2(:,2),partibib2(:,3),partibib2(:,4))
title('Max Slope Method, Bare Wire Embedded Thermocouple, Ice Water to Boiling Water')
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-1 1 -inf 1])
subplot(3,1,2)
plot(partiibib2(:,1),partiibib2(:,2),partiibib2(:,3),partiibib2(:,4))
xlabel('Time (s)')
ylabel('Temperature (�C)')
legend('Data','Prediction','location','northeast')
axis([-1 1 -10 110])
grid on
subplot(3,1,3)
plot(partiiibib2(:,1),partiiibib2(:,2))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
legend ('Data Minus Prediction','location','southeast')
axis([-1 1 -20 20])

%% Dynamic Calibration Part 4

%finding Syx values
sbisyx = Syx(steelboilicearray(:,1),steelboilicearray(:,2));
sbiss = num2str(sbisyx,3);
sbist = strcat('Syx = ',sbiss);

sbisyx2 = Syx(steelboilicearray2(:,1),steelboilicearray2(:,2));
sbiss2 = num2str(sbisyx2,3);
sbist2 = strcat('Syx = ',sbiss2);

sibsyx = Syx(steeliceboilarray(:,1),steeliceboilarray(:,2));
sibss = num2str(sibsyx,3);
sibst = strcat('Syx = ',sibss);

sibsyx2 = Syx(steeliceboilarray(:,1),steeliceboilarray(:,2));
sibss2 = num2str(sibsyx2,3);
sibst2 = strcat('Syx = ',sibss2);

abisyx = Syx(alumboilicearray(:,1),alumboilicearray(:,2));
abiss = num2str(abisyx,3);
abist = strcat('Syx = ',abiss);

abisyx2 = Syx(alumboilicearray2(:,1),alumboilicearray2(:,2));
abiss2 = num2str(abisyx2,3);
abist2 = strcat('Syx = ',abiss2);

aibsyx = Syx(alumiceboilarray(:,1),alumiceboilarray(:,2));
aibss = num2str(aibsyx,3);
aibst = strcat('Syx = ',aibss);

aibsyx2 = Syx(alumiceboilarray2(:,1),alumiceboilarray2(:,2));
aibss2 = num2str(aibsyx2,3);
aibst2 = strcat('Syx = ',aibss2);

bbisyx = Syx(bareboilicearray(:,1),bareboilicearray(:,2));
bbiss = num2str(bbisyx,3);
bbist = strcat('Syx = ',bbiss);

bbisyx2 = Syx(bareboilicearray2(:,1),bareboilicearray2(:,2));
bbiss2 = num2str(bbisyx2,3);
bbist2 = strcat('Syx = ',bbiss2);

bibsyx = Syx(bareiceboilarray(:,1),bareiceboilarray(:,2));
bibss = num2str(bibsyx,3);
bibst = strcat('Syx = ',bibss);

bibsyx2 = Syx(bareiceboilarray2(:,1),bareiceboilarray2(:,2));
bibss2 = num2str(bibsyx2,3);
bibst2 = strcat('Syx = ',bibss2);

%plot residuals on top of each other
%steel boil ice residuals
figure(41)
subplot(3,1,1)
plot(D3sbi(:,5),D3sbi(:,6))
title('Residuals - Boiling Water to Ice Water, 5\sigma Method')
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,sbist)
legend ('Steel Embedded Thermocouple','location','southeast')
%aluminum boil ice residuals
subplot(3,1,2)
plot(D3abi(:,5),D3abi(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,abist)
legend ('Aluminum Embedded Thermocouple','location','southeast')
%bare boil ice residuals
subplot(3,1,3)
plot(D3bbi(:,5),D3bbi(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,bbist)
legend ('Bare Wire Thermocouple','location','southeast')

figure(42)
%steel ice boil residuals
subplot(3,1,1)
plot(D3sib(:,5),D3sib(:,6))
title('Residuals - Ice Water to Boiling Water, 5\sigma Method')
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,sibst)
legend ('Steel Embedded Thermocouple','location','southeast')

%aluminum ice boil residuals
subplot(3,1,2)
plot(D3aib(:,5),D3aib(:,5))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,aibst)
legend ('Aluminum Embedded Thermocouple','location','southeast')

%bare ice boil residuals
subplot(3,1,3)
plot(D3bib(:,5),D3bib(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,bibst)
legend ('Bare Wire Embedded Thermocouple','location','southeast')

figure(43) %%%%%%%%%%%%%%%%%max slope method
subplot(3,1,1)
plot(D3sbi2(:,5),D3sbi2(:,6))
title('Residuals - Boiling Water to Ice Water, Max Slope Method')
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,sbist2)
legend ('Steel Embedded Thermocouple','location','southeast')
%aluminum boil ice residuals
subplot(3,1,2)
plot(D3abi2(:,5),D3abi2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,abist2)
legend ('Aluminum Embedded Thermocouple','location','southeast')
%bare boil ice residuals
subplot(3,1,3)
plot(D3bbi2(:,5),D3bbi2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,bbist)
legend ('Bare Wire Embedded Thermocouple','location','southeast')

figure(44)
%steel ice boil residuals
subplot(3,1,1)
plot(D3sib2(:,5),D3sib2(:,6))
title('Residuals - Ice Water to Boiling Water, Max Slope Method')
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,sibst2)
legend ('Steel Embedded Thermocouple','location','southeast')

%aluminum ice boil residuals
subplot(3,1,2)
plot(D3aib2(:,5),D3aib2(:,5))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,aibst2)
legend ('Aluminum Embedded Thermocouple','location','southeast')

%bare ice boil residuals
subplot(3,1,3)
plot(D3bib2(:,5),D3bib2(:,6))
xlabel('Time (s)')
ylabel('Residuals (�C)')
grid on
axis([-5 30 -20 20])
text(-4.5,10,bibst2)
legend ('Bare Wire Embedded Thermocouple','location','southeast')







%% Dynamic Calibration Part 5
load lab2part1variables.mat

bareiceairtime = xlsread('Michalak_Popecki_Rose.xlsx',5,'a:a'); %this is for bareiceair3 - three samples were taken and this one has the best data
bareiceairvoltage = xlsread('Michalak_Popecki_Rose.xlsx',5,'B9:B12008');
%a very noisy signal...

bareiceairvoltage = smooth(bareiceairvoltage,1001);
bareiceairarray = pros(bareiceairtime,bareiceairvoltage,2.19);
%bareiceairtemperature = ((bareiceairarray(:,2))-betaHat(1))/betaHat(2); %betahat 2 is the slope

rt = 26; %�C, from my lab notebook
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

%% Funtions used

% function [output] = middlefit(xdata,ydata) %x input is time, y input is temperature from data
% the inputs are shifted such that the 0 value of the xdata is the beginning
% of the event.
% 
% this function provides an output array containing the information to be
% plotted on the middle plot of the 3-plot subplot
% 
% gamma=(ydata(end)-ydata)/(ydata(end)-ydata(1));
% for i=1:length(xdata)
%     if gamma(i)<0.05
%         endgamma=i;
%         break
%     end
% end
% 
% lngamma=log(gamma(1:endgamma));
% num=lngamma.*xdata(1:endgamma);
% den=xdata(1:endgamma).^2;
% ao=sum(num)/sum(den);
% 
% 
% Tfinal = ydata(end); %the final temperature of the data
% 
% Tinitloc = find(xdata==0); %the location in the array where time is zero
% 
% Tinitial = ydata(Tinitloc); %the Temperature value where time is zero
% 
% tau = (1/ao)*-1;
% 
% Tpred = Tfinal - (Tfinal-Tinitial)*exp(-xdata/tau);
% 
% 
% output = [xdata,ydata,xdata,Tpred]; %[data time, data Temperature, predicted time (same as data time), predicted Temperature]
% 
% 
% end
% 
% function [output] = Syx(xdata,ydata)
% Computers the Syx of a data set
% 
% tstart = find(xdata==0,1);
% xnew = xdata(tstart:length(xdata));
% ynew = ydata(tstart:length(ydata));
% 
% fitdata = polyfit(xnew,ynew,1);
% yfromline = fitdata(1)*xnew+fitdata(2);
% 
% Yc = yfromline; %The value of y predicted by the polynomial equation for a given value of x
% tvp = 2.262; %For N = 10, 95% confidence
% yiyci = (ynew-yfromline).^2;
% sumyiyci = sum(yiyci);
% Syx = (sumyiyci/(length(ynew)-1))^.5; %standard error of the fit
% 
% disp(Syx)
% output = Syx;
% 
% end
% 
% function [output] = p2(xdata,ydata)
% this function outputs an array with ready to plot info for the two curves
% and the residuals
% 
% where xdata is the time and ydata is the temperature
% 
% Tfinal = ydata(end); %the final temperature of the data
% Tinitloc = find(xdata==0); %the location in the array where time is zero
% Tinitial = ydata(Tinitloc); %the Temperature value where time is zero
% Tattau = Tinitial+.632*(Tfinal - Tinitial); %this is the temperature at the point where t is equal to tau
% disp(Tattau)
% 
% find the time at which the temperature is equal to Tattau - this is our
% time constant tau
% 
% conditional statement 
% if Tfinal<Tinitial
%     tauindex = find(ydata<=Tattau,1); %the location of tau in the array
% else
%     tauindex = find(ydata>=Tattau,1);
% end
%     
% 
% tau = xdata(tauindex);
% 
% Tpred = Tfinal-(Tfinal-Tinitial)*exp(-xdata/tau);
% 
% residuals = ydata-Tpred;
% 
% output = [xdata,ydata,xdata,Tpred,xdata,residuals]; %[data time, data Temperature, predicted time (same as data time), predicted Temperature]
% 
% end
% 
% function [output] = pros(time,voltage,startingtime)
% load lab2part1variables.mat
% offset = 0;
% 
% for i=1:length(time)
%     if time(i)>startingtime
%         basetime=i;
%         break
%     end
% end
% baseline=mean(voltage(1:basetime)); %average of data in baseline region
% basedev=std(voltage(1:basetime)); %standard deviation of baseline region
% threshold=5*basedev; %threshold to define the start of an event - 5 sigma
% 
% for i=1:length(time)
%     if (abs(voltage(i)-baseline)>threshold)
%         starttime=i;
%         break
%     end
% end
% create new variables that start from t = 0 and only contain event data
% newtime=time(starttime:length(time))-time(starttime); %the commented out portion starts time t=0 at the event initiation point
% time = time-time(starttime); %shifts time such that t=0 
% newvoltage=voltage(starttime:length(time));
% 
% Tstart = time(starttime);
% 
% Vstart = voltage(starttime);
% 
% when it is more useful to output a temperature instead of voltage
% tcv = (voltage*1000-betaHat(1))/betaHat(2); %�C newvoltage for trimmed, voltage for untrimmed
% 
% 
% output = [tout,vout];
% output = [time,tcv]; %trimmed: newtime,newvoltage. untrimmed: time,voltage
% end
% 
% function [output] = slide(time,voltage)
% outputs the maximum slope and its position in the matrix
% time is x, voltage is y
% load lab2part1variables.mat
% 
% for i = 1:1:(length(time)-51)
%     limit = (i+51);
%     tmask = time(i:limit);
%     vmask = voltage(i:limit);
%     fit = polyfit(tmask,vmask,1);
%     m = fit(1);
%     slopes(i) = m; %an array of slopes at each point in the line
%     posslopes = abs(slopes);
%     [~,pos] = max(posslopes);
%     maxslope = slopes(pos);
%     
%     newtime = time(pos:length(time));
%     time = time-time(pos);
%     newvoltage = voltage(pos:length(voltage));
%     
%     tcv = (voltage*1000-betaHat(1))/betaHat(2); %�C newvoltage for trimmed, voltage for untrimmed
%     
%     output = [time,tcv];
% end
% 
% 
% 
% end
% 
% function [output] = gammafit(xdata,ydata) %x is time
% gamma=(ydata(end)-ydata)/(ydata(end)-ydata(1));
% for i=1:length(xdata)
%     if gamma(i)<0.05
%         endgamma=i;
%         break
%     end
% end
% 
% lngamma=log(gamma(1:endgamma));
% num=lngamma.*xdata(1:endgamma);
% den=xdata(1:endgamma).^2;
% ao=sum(num)/sum(den);
% predictln=ao*xdata(1:endgamma);
% figure(4000)
% plot(xdata(1:endgamma),lngamma,xdata(1:endgamma),predictln)
% 
% output = [xdata(1:endgamma),lngamma,xdata(1:endgamma),predictln];
% 
% 
% 
% end
% 
% function [output] = bottomfit(middlefit)
% This provides an array to plot the residuals vs. time
% 
% the residuals are the difference between the data and the prediction
% this function will use the results of "middlefit.m" as an input
% 
% datay = middlefit(:,2); %the second column of the middlefit output array
% predy = middlefit(:,4); %the predicted y values
% xdata = middlefit(:,1); %remember that the xvalues are the same for both arrays
% 
% residuals = datay - predy; %an array of differences between the actual data and the prediction
% 
% output = [xdata,residuals];
% 
% 
% 
% end
% 
























