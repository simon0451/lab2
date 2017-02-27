load lab2variables.mat

figure(1)
plot(partiabi2(:,1),partiabi2(:,2),partiabi2(:,3),partiabi2(:,4))
xlabel('Time (s)')
ylabel('ln(\Gamma)')
legend('Data','Prediction','location','northeast')
grid on
axis([-5 20 -inf 1])