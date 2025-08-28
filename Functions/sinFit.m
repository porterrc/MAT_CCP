function [fitresult, gof] = sinFit(cellbaz, cellamp)
%CREATEFIT(CELLBAZ,CELLAMP)

%% Fit: 'untitled fit 1'.
cellbaz = cellbaz(~isnan(cellamp));
cellbaz = cellbaz*pi/180;
cellamp = cellamp(~isnan(cellamp));
cellamp = detrend(cellamp,'constant');


[xData, yData] = prepareCurveData( cellbaz,cellamp );

% Set up fittype and options.
ft = fittype( 'a*sin(x+b)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -0];
opts.StartPoint = [0 0];
opts.Upper = [2 2*pi];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );


