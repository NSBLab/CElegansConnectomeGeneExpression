function [h, yMeans, xThresholds]  = BF_PlotQuantiles(xData,yData,numThresholds,alsoScatter,makeNewFigure,theColor,addStd)
% Plots x-y scatter, but with mean of y plotted in quantiles of x
% Ben Fulcher
%-------------------------------------------------------------------------------

if nargin < 3 || isempty(numThresholds)
    numThresholds = 10;
end
if nargin < 4
    alsoScatter = 0;
end
if nargin < 5
    makeNewFigure = false;
end
if nargin < 6
    theColor = 'k';
end
if nargin < 7
    addStd = false;
end

%-------------------------------------------------------------------------------
% Filter out NaNs:
goodBoth = (~isnan(xData) & ~isnan(yData));
if ~any(goodBoth)
    error('No good data');
elseif any(~goodBoth)
    xData = xData(goodBoth);
    yData = yData(goodBoth);
    fprintf(1,'Removed %u bad samples from x/y data\n',sum(~goodBoth));
end

xThresholds = arrayfun(@(x)quantile(xData,x),linspace(0,1,numThresholds));
xThresholds(end) = xThresholds(end) + eps; % make sure all data included in final bin
yMeans = arrayfun(@(x)mean(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
if addStd
    yStds = arrayfun(@(x)std(yData(xData>=xThresholds(x) & xData < xThresholds(x+1))),1:numThresholds-1);
end

% ------------------------------------------------------------------------------
% Plot:
if makeNewFigure
    f = figure('color','w'); box('on'); hold on
end
theStyle = '-';
lineWidthMajor = 2;
lineWidthMinor = 1;

if alsoScatter
    plot(xData,yData,'.','color',brighten(theColor,0.4));
end

for k = 1:numThresholds-1
    h = plot(xThresholds(k:k+1),ones(2,1)*yMeans(k),'LineStyle',theStyle,'LineWidth',lineWidthMajor,'Color',theColor);
    if addStd
        plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)+yStds(k)),'LineStyle',':','LineWidth',lineWidthMinor,'Color',theColor)
        plot(xThresholds(k:k+1),ones(2,1)*(yMeans(k)-yStds(k)),'LineStyle',':','LineWidth',lineWidthMinor,'Color',theColor)
        if k > 1
            % Also add joiner from previous error bar:
            plot(xThresholds(k)*ones(2,1),[yMeans(k-1)+yStds(k-1),yMeans(k)+yStds(k)],'LineStyle',':','LineWidth',lineWidthMinor,'Color',theColor)
            plot(xThresholds(k)*ones(2,1),[yMeans(k-1)-yStds(k-1),yMeans(k)-yStds(k)],'LineStyle',':','LineWidth',lineWidthMinor,'Color',theColor)
        end
    end
    plot(mean(xThresholds(k:k+1)),yMeans(k),'o','MarkerSize',5,'LineStyle',theStyle,'LineWidth',lineWidthMajor,'Color',theColor)
end

end
