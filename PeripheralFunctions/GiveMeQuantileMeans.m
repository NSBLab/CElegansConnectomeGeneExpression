% Outputs the xRange used for the quantiles
% Outputs the means of the yData provided in each quantiled-bin
% Outputs the indices of the data used in each quantile (as a cell)
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-11-05


function [xRange,yMeans,affectedInd] = GiveMeQuantileMeans(xData,yData,numThresholds)

% Set up quantiles for computing distances:
xQuantiles = linspace(0,1,numThresholds + 1);

% xRange = zeros(numThresholds+1,1);
for i = 1:numThresholds+1
    xRange(i) = quantile(xData,xQuantiles(i));
end
xRange(1) = xRange(1) - eps;

% affectedInd contains the indices for links in the given distance range (by quantile)
affectedInd = cell(numThresholds,1);
yMeans = zeros(numThresholds,1);
for i = 1:numThresholds
    affectedInd{i} = (xData > xRange(i) & xData <= xRange(i+1));
    if ~isempty(yData)
        yMeans(i) = nanmean(yData(affectedInd{i}));
    end
end

    
end
