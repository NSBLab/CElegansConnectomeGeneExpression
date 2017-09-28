function hubConnectionDistance(C,plotWhat)
%===============================================================================
%===============================================================================
% Idea is to compare connections involving hubs to non-hubs in each anatomical type
% cf probabilityVSdistance
%===============================================================================
%===============================================================================

if nargin < 2
    plotWhat = {'violin','bar'};
end

%-------------------------------------------------------------------------------
Adj = GiveMeAdj(C,'zeroBinary');
[~,~, deg] = degrees_dir(Adj);
dist = GiveMeDist(C);

%-------------------------------------------------------------------------------
maskAnatomy = GiveMeMask(C,'AnatomyPairs');
theMasks = fieldnames(maskAnatomy);

D = GiveMeDefault();
isHub = double(deg > D.kHub);
isRich = isHub'*isHub;
isRich(logical(eye(size(isRich)))) = 0;
notRich = double(~isRich);
notRich(logical(eye(size(notRich)))) = 0;

%===============================================================================
if ismember('violin',plotWhat);
    f = figure('color','w');
    for k = 1:length(theMasks)
        ax = subplot(3,3,k);
        maskedD = (maskAnatomy.(theMasks{k}) & Adj);
        dRich = dist(maskedD & isRich);
        dNotRich = dist(maskedD & notRich);
        JitteredParallelScatter({dRich,dNotRich},true,true,false);
        ax.XTick = [1,2];
        ax.XTickLabel = {'rich','notRich'};
        ylabel('Distance (mm)')

        if ~isempty(dRich)
            [~,p] = ttest2(dRich,dNotRich,'Vartype','unequal');
            title(sprintf('%s: p=%.3g',theMasks{k},p))
        else
            title(theMasks{k})
        end
    end
end

%-------------------------------------------------------------------------------
% But does the proportion change?
if ismember('bar',plotWhat);
    numConns = zeros(length(theMasks),1);
    for k = 1:length(theMasks)
        hubConns = maskAnatomy.(theMasks{k}) & Adj & isRich;
        nonhubConns = maskAnatomy.(theMasks{k}) & Adj & notRich;
        numConns(k,1) = sum(hubConns(:));
        numConns(k,2) = sum(nonhubConns(:));
    end
    numConnsNorm = numConns;
    numConnsNorm(:,1) = numConnsNorm(:,1)/sum(numConnsNorm(:,1));
    numConnsNorm(:,2) = numConnsNorm(:,2)/sum(numConnsNorm(:,2));

    f = figure('color','w'); ax = gca;
    bar(numConnsNorm,'grouped','Horizontal','on')
    legend({'rich','non-rich'})
    ax.YTick = 1:length(theMasks);
    theMasksHyphen = cellfun(@(x)[x(1:4),'->',x(5:end)],theMasks,'UniformOutput',false);
    ax.YTickLabel = theMasksHyphen;
    [rgb_colorMatrix,labels] = GiveMeColors('richFeederPeripheral');
    colormap([rgb_colorMatrix(4,:);ones(1,3)*0.5]);
    xlabel(ax,'Proportion of connections')
    f.Position = [1479         980         417         249];
end
