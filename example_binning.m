
clear
close all
clc

%% generate example noisy data
npts=2000;
xdata=linspace(0,2*pi,npts);
ydata=cos(xdata) + 2*randn(1,npts);


%% bin data with bootstrapped confidence intervals
func=@nanmean;
Nbins=25;

[xbin,ybin]=mydiscretize_CI(xdata,ydata,func,Nbins);


%% plot

figure
scatter(xdata,ydata,'k.')
hold all
errorbar_CI(xbin,ybin,'ro')
plot(xdata,cos(xdata),'b-')

legend('data with noise','binned noisy data','true signal')