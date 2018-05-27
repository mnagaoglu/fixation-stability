% demo script for computing fixation stability using ISOA method
%
% MNA 5/26/2018 wrote it. mnagaoglu@gmail.com
%
close all;
clc;


isShowPlot = 1; % to see the figures
cumulativeProbability = 0.68; % must be between 0 and 1.

% number of data points
simN = 10000;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple example first. Normally distributed data with some correlation
xDeg = (1 + .5*randn(simN,1));
yDeg = (-2 + 1*randn(simN,1)) + 0.3*xDeg;

% compute density, PRL, ISOA, and even BCEA for comparison
[isoa, bcea, PRL, PRL2, density, xGrid, yGrid, fh] = ...
    ComputeFixationStability(xDeg, yDeg, cumulativeProbability, isShowPlot); %#ok<*ASGLU>
fh.Name = 'Example I';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% non-Gaussian distribution (which might occur if subject fixates 
% different locations creating multiple "islands" of data points)
xDeg = [(1 + .5*randn(simN/2,1)); (-2+ 2*randn(simN/2,1))];
yDeg = [(-2 + 1*randn(simN/2,1)); (2+ 1*randn(simN/2,1))] + 0.3*xDeg;

% compute density, PRL, ISOA, and even BCEA for comparison
[isoa, bcea, PRL, PRL2, density, xGrid, yGrid, fh] = ...
    ComputeFixationStability(xDeg, yDeg, cumulativeProbability, isShowPlot);
fh.Name = 'Example II';