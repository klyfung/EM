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

% Extracting the edge.
% Define the start of the edge.
startedge = xdata > 176;
xdata1 = xdata(startedge);
ydata1 = ydata(startedge);
% Define end of the edge.
endedge = xdata1 < 381;
xdata2 = xdata1(endedge);
ydata2 = ydata1(endedge);

%% Define colour order for residuals plots
% Note, second entry in matrix is first colour of plot.
co = [0.00 0.00 0.00
	1.00  0.40  0.40
	1.00  0.70  0.40
	1.00  1.00  0.40
	0.70  1.00  0.40
	0.40  1.00  0.40
	0.40  1.00  0.70
	0.40  1.00  1.00
	0.40  0.70  1.00
	0.40  0.40  1.00
	0.70  0.40  1.00
	1.00  0.40  1.00
	1.00  0.40  0.70];
set(groot,'defaultAxesColorOrder',co)

%% Fitting 'for' loop
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

% (i = start:increment:end for excluded data points)
for i = 200:10:280
exclude1 = xdata2 > i;
% Define fit used and 'Exclude' data points in xdata (eV) for fitting
f = fit(xdata2,ydata2,'power2','Exclude',exclude1);
% Get residuals from fit to plot subtracted spectra
residuals = ydata2 - f(xdata2);

% Plot subtracted spectra
p1 = plot(xdata2,residuals);
ax1 = gca;
% Define colours of plots and names for legend
colorOrder = get(gca, 'ColorOrder');
set(p1,'Color',colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :),'LineWidth',2,'DisplayName',num2str(i))
hold on
end

% Plot original EELS data
plot(xdata2,ydata2,'Color','k','LineWidth',2,'DisplayName','Original EELS data');
hold off

%% Define characteristics for axes
ax1.XLim = [-inf inf]; % Limits of x-axis
ax1.YLim = [-inf inf]; % Limits of y-axis
ax1.FontName = 'Calibri';
ax1.FontSize = 30;
ax1.TickDir = 'out';
ax1.TickLength = [0.005 0.005];
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.Layer = 'bottom';
ax1.Title.String = 'EEL spectra after subtracting fitted curves';
ax1.Title.FontWeight = 'normal';
ax1.XLabel.String = 'eV';
ax1.YLabel.String = 'Counts';
lgd1 = legend(ax1,{},'FontSize',30,'FontWeight','normal','box','off','Location','Northeastoutside');
title(lgd1,'Fits excluding data above (i) eV','FontSize',30,'FontWeight','normal')