%% Plotting background subtracted EEL spectra
% KLYF 2020
% Published in https://www.sciencedirect.com/science/article/pii/S0304399120302035

clc
close all
clear all

% To extract and plot the fitted EELS spectrum, save xdata2 (eV) and
% residuals (counts after background subtraction) from the Workspace after
% running EELS_fit_analysis.m using a value of (i) that gives the best fit.

% Import original .msa data exported from Digital Micrograph 3.
filename = 'file_name.msa'; % This is your file name.
delimiterIn = ','; % This is the character that separates the two columns of data.
headerlinesIn = 20; % This is the number of lines of text at the start of the data that are skipped.
msadata = importdata(filename,delimiterIn,headerlinesIn);
data = msadata.data;

load ('EELS_fit_xdata_file_name.mat')
load ('EELS_fit_ydata_residuals_file_name.mat')

% Assign variables from imported data.
x1 = data(:,1);
y1 = data(:,2);
x2 = xdata2;
y2 = residuals;

% Saves subtracted spectrum workspace variables from MATLAB to .txt format
% for plotting in other programs.
t1 = table(x2,y2);
writetable(t1,'subtracted-spectrum.txt','WriteRowNames',true)

% Plot spectrum
p1 = plot(x1,y1);
ax1 = gca;
set(p1,'Color','k','LineWidth',2,'DisplayName','Original EEL spectrum')
hold on
p2 = plot(x2,y2);
set(p2,'Color','r','LineWidth',2,'DisplayName','EEL spectrum after fitting')
hold off

%% Define characteristics of axes
ax1.XLim = [-inf inf];
ax1.YLim = [-inf inf];
ax1.FontName = 'Calibri';
ax1.FontSize = 30;
ax1.TickDir = 'out';
ax1.TickLength = [0.005 0.005];
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.Layer = 'bottom';
ax1.Title.String = 'Original and subtracted EEL spectra';
ax1.Title.FontWeight = 'normal';
ax1.XLabel.String = 'eV';
ax1.YLabel.String = 'Counts';
lgd1 = legend(ax1,{},'FontSize',30,'FontWeight','normal','box','off','Location','Northeast');
title(lgd1,[],'FontSize',30,'FontWeight','normal')