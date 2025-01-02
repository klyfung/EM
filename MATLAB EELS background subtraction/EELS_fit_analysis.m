%% EELS background nonlinear least squares fitting (using the Curve Fitting Toolbox)
% KLYF 2020
% Published in https://www.sciencedirect.com/science/article/pii/S0304399120302035

clc
close all
clear all

% Import .msa data exported from Digital Micrograph 3.
filename = 'file_name.msa'; % This is your file name.
delimiterIn = ','; % This is the character that separates the two columns of data.
headerlinesIn = 20; % This is the number of lines of text at the start of the data that are skipped.
msadata = importdata(filename,delimiterIn,headerlinesIn);
data = msadata.data;

% Assign variables from imported data (xdata = ev; ydata = counts).
xdata = data(:,1);
ydata = data(:,2);

% There may be situations where the data of interest is a smaller section
% of the total spectrum. The code below allows the user to define the start
% and end of that data. A larger background before the start of the edge
% tends to lead to better fitting of the background.

% Extracting the edge. Define the start of the edge.
startedge = xdata > 176;
xdata1 = xdata(startedge);
ydata1 = ydata(startedge);
% Define end of the edge.
endedge = xdata1 < 381;
xdata2 = xdata1(endedge);
ydata2 = ydata1(endedge);

% Create table of imported data for exclusion of data for fitting.
table1 = table(xdata2,ydata2);

%% Analysing exponential fitting
% Exclude data from table where xdata is above i eV and change i to the
% appropriate value.
i = 280;
data2 = table1{table1.xdata2 < i,:};
% Make new variables from excluded data for fitting.
xdata3 = data2(:,1);
ydata3 = data2(:,2);

% Define fit - use one of the fits in single quotations below:

% 'exp1' is a one-term exponential model f(x) = a*exp(b*x)

% 'exp2' is a two-term exponential model f(x) = a*exp(b*x) + c*exp(d*x)

% 'power1' is a one-term power model f(x) = a*x^b

% 'power2' is a two-term power model f(x) = a*x^b + c

% If any other fits are needed, you need to define the fittype as they're
% not included in the MATLAB Curve Fitting Toolbox e.g. exp3 is a
% three-term exponential model f(x) = a*exp(b*x) + c*exp(d*x) + e*exp(f*x).
% exp3 = fittype('a*exp(b*x) + c*exp(d*x)+
% e*exp(f*x)','dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f'});
% You then need to define reasonable start points if exp3 is to be used and
% input exp3 and not 'exp3'. This type of fitting wasn't done as part of
% this script because defining the starting points of six different
% coefficients gets rather complicated and, for the most part, the fits
% available in the Curve Fitting Toolbox are sufficient.

% Fit the new excluded data xdata2 (eV) - the fit in single quotations can
% be changed as appropriate.
[f,gof,output] = fit(xdata3,ydata3,'power2')
% Get residuals from fit
residuals = ydata2 - f(xdata2);
% Get fit options
fitoptions = fitoptions(f)

% Approximate the signal integral using the trapezoidal rule. This can be
% used to quantify the signal-to-noise ratio (SNR) as well as in further
% quantification of chemical composition. The partial inelastic scattering
% cross-section needed for absolute quantification depends on the system
% being studied. The window of integration can be defined by the user.
% Define the start of integration for signal integral.
startint = xdata2 > 284;
xdataint1 = xdata2(startint);
residualsint1 = residuals(startint);
% Define end of integration for signal integral.
endint = xdataint1 < 300;
xdataint2 = xdataint1(endint);
residualsint2 = residualsint1(endint);
% Integrate signal
ik = trapz(xdataint2,residualsint2)
% Integrate background
fityvalues = f(xdata2);
bkgint1 = fityvalues(startint);
bkgint2 = bkgint1(endint);
ib = trapz(xdataint2,bkgint2)
% Calculate variance in the background integral
varib = var(bkgint2)
% h parameter
h = (ib+varib)/ib
% Signal-to-noise ratio (SNR)
snr = ik/((ik+(h*ib))^0.5)

% Plot fit with confidence bounds, original EELS data, and the subtracted
% spectrum. The confidence bounds for the fitted coefficients determine the
% accuracy of the fit. Bounds that are far apart indicate uncertainty in
% the fit.
subplot(1,2,1)
% Plot fit with confidence bounds and original data
p1 = plot(f,'b',xdata2,ydata2,'k','predobs');
hold on
% Plot subtracted spectrum
p2 = plot(xdata2,residuals,'r');
hold on
ax1 = gca;
% Define characteristics of plot
set(p1,'LineWidth',2)
set(p2,'LineWidth',2)

% Residuals from a fitted model are the differences between the original
% data and the fit to the original data. Assuming the model is correct, the
% residuals approximate the random errors. The model fits the data well if
% the residuals appear to behave randomly. If the residuals have a
% systematic pattern, the model does not fit the data very well. Plot
% residuals
subplot(1,2,2)
p3 = plot(f,xdata3,ydata3,'Residuals');
hold off
ax2 = gca;
% Define characteristics of plot
set(p3,'LineWidth',2)

%% Define characteristics of fit and original data axes
ax1.XLim = [-inf inf]; % Limits of x-axis
ax1.YLim = [-inf 500000]; % Limits of y-axis
ax1.FontName = 'Calibri';
ax1.FontSize = 30;
ax1.TickDir = 'out';
ax1.TickLength = [0.005 0.005];
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.Layer = 'bottom';
ax1.Title.String = '';
ax1.Title.FontWeight = 'normal';
ax1.XLabel.String = 'eV';
ax1.YLabel.String = 'Counts';
lgd1 = legend(ax1,{'Original data','Fitted curve','Upper prediction bounds','Lower prediction bounds','Subtracted spectrum'},'FontSize',30,'FontWeight','normal','box','off','Location','North');
title(lgd1,['Background fit of data below ' num2str(i) ' eV'],'FontSize',30,'FontWeight','normal')

%% Define characteristics for residuals axes
ax2.XLim = [-inf inf]; % Limits of x-axis
ax2.YLim = [-inf inf]; % Limits of y-axis
ax2.FontName = 'Calibri';
ax2.FontSize = 30;
ax2.TickDir = 'out';
ax2.TickLength = [0.005 0.005];
ax2.XGrid = 'on';
ax2.YGrid = 'on';
ax2.Layer = 'bottom';
ax2.Title.String = '';
ax2.Title.FontWeight = 'normal';
ax2.XLabel.String = 'eV';
ax2.YLabel.String = 'Counts';
lgd2 = legend(ax2,{'Residuals','Zero line'},'FontSize',30,'FontWeight','normal','box','off','Location','North');
title(lgd2,'Residuals of fit','FontSize',30,'FontWeight','normal')