function [S,P] = RichClub(C,G,analyzeWhat,measure,coexpMeasure,whatMask,whatAdj,networkType,useOnlyInterneurons,doRanksum)
% ------------------------------------------------------------------------------
% Function plots coexpression for rich/feeder/peripheral lins as a function
% of degree using mean to summarise coexpression at each threshold
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% G (gene data structure)
% analyzeWhat 'coexpression', 'lineage distance', 'connection distance';
% coexpMeasure - choose coexpression measure as defined in G.Corr. default Pearson_noLR
% ------------------------------------------------------------------------------

if nargin < 3
    analyzeWhat = 'coexpression';
end
if nargin < 4 || isempty(measure)
    measure = 'median';
end
if nargin < 5 && strcmp(analyzeWhat, 'coexpression')
    coexpMeasure = 'Pearson';
    fprintf(1,'Using Pearson without LR coexpression BY DEFAULT\n');
end
if nargin == 5 && strcmp(analyzeWhat,'coexpression')
    coexpMeasure=coexpMeasure;
end
if nargin < 6
    whatMask = 'connected'; % 'connected', 'connectedEither', 'rich'
end
if nargin < 7
    whatAdj = 'zeroBinary'; % 'rawBinary','raw weighted','zero binary','zero weighted'
end
if nargin < 8
    networkType = 'bd'; % 'bu','bd','wu','wd'
end
if nargin < 9
    useOnlyInterneurons = false; % 1-use only interneurons; 2 - use all neurons except interneurons or 0-use all neurons
end
if nargin < 10
    doRanksum = true;
end
pThreshold = 0.05;
D = GiveMeDefault();
% ------------------------------------------------------------------------------
%% INPUTS:
% ------------------------------------------------------------------------------

% ------------- (1) analyzeWhat: what to plot the rich club curves for:

switch analyzeWhat
    case 'coexpression'
        realLinkData = GiveMeCoexpression(G,coexpMeasure,true);
    case 'connection distance'
        realLinkData = GiveMeDist(C,[],true);
    case 'lineage distance'
        realLinkData = C.LineageDistance_noLR;
    case 'birth time difference'
        realLinkData = C.BirthTimeDiff_noLR;
        %         bornEarly = C.BirthTime<1000;
        %         earlyMask = bornEarly.*bornEarly';
        %         realLinkData = C.BirthTimeDiff_noLR.*earlyMask; %C.BirthTimeDiff_noLR;
        %realLinkData(realLinkData>700)=nan;
end

% -------------- (3) labelNodesHow: how to define "rich"
labelNodesHow = 'hub-kth'; extraParam = {'degree',0};

numBins = 'all'; % Range of k to plot the rich club curve across

plotDist = false; % plot full distributions (at each k)

% ------------------------------------------------------------------------------
%% Assign data measured at each link in the network
% ------------------------------------------------------------------------------
% All directed connections:
linkedAdj = GiveMeAdj(C,whatAdj);
% Add a mask to only include particular types of connections
mask = GiveMeMask(C, whatMask,networkType,linkedAdj,extraParam,1);
numNeurons = size(linkedAdj,1);
linkedAdj = mask.special;

nn = linspace(1,numNeurons,numNeurons);
inter = find(C.RegionM(:,10)==1);
ninter = setdiff(nn, inter);
if useOnlyInterneurons == 1
    realAdj.linked = linkedAdj(inter, inter);
    allLinkData = realLinkData(inter, inter);
elseif useOnlyInterneurons == 2
    realAdj.linked = linkedAdj(ninter, ninter);
    allLinkData = realLinkData(ninter, ninter);
else
    realAdj.linked = linkedAdj;
    allLinkData = realLinkData;
end

allLinkData(realAdj.linked==0) = NaN;

numNodes = length(realAdj.linked);

% ------------------------------------------------------------------------------
% nodeData (~degree) should not change with different nulls, which preserve the in/out degree of all nodes
% ------------------------------------------------------------------------------
[~,~,nodeData] = AdjLabelNodes(labelNodesHow,realAdj.linked,extraParam,networkType);

% ------------------------------------------------------------------------------
% Get groups of links based on their degree (or use all k)
% ------------------------------------------------------------------------------
sortK = sort(nodeData,'descend');
maxK = sortK(2); % Up to the second-highest k
if strcmp(numBins,'all')
    % All k are a bin
    kr = min(nodeData):maxK;
else % numeric number of bins
    kr = linspace(min(nodeData),maxK,numBins);
    kr = round(kr);
end
krAll = min(nodeData):max(nodeData);

% Proportion of nodes called a 'hub' through kr
switch measure
    case 'mean'
        propisHub = arrayfun(@(x) mean(nodeData > x),kr);
    case 'median'
        propisHub = arrayfun(@(x) median(nodeData > x),kr);
end

% ------------------------------------------------------------------------------
% Go through each class of links and compute statistics on the set of link
% data compared to the nulls:

whatLinks = {'rich','feedin','feedout','local'};
%whatLinks = {'rich','feedboth','local'};

% Ben Fulcher, 2015-01-06
% Computes a t-test between special and non-special links for each k, and each link-type:
tStats = zeros(length(kr),length(whatLinks));
allHubHub = cell(1,1); % make a 1-component cell for consistency with null version
allHubHub{1} = cell(length(kr),length(whatLinks));
for i = 1:length(kr)
    for j = 1:length(whatLinks) % loop across rich, feedin, feedout, and local connections:

        % Keep links between high degree nodes that have linkData (i.e., links that exist):
        % (assumption that nodeData will be the same for all null networks)
        % Ben Fulcher, 2015-01-06 -- changed from >= to > to match actual rich-club coeff definition
        r = (nodeData > kr(i));

        keepMe = logical(zeros(numNodes,numNodes));
        % Add ones for a given type of link:
        switch j
            case 1 % 'rich'
                keepMe(r,r) = 1;
            case 2 % 'feedin'
                keepMe(~r,r) = 1;
            case 3 % 'feedout'
                keepMe(r,~r) = 1;
                %             case 4 % 'feedboth'
                %                 keepMe(r,~r) = 1;
                %                 keepMe(~r,r) = 1;
            case 4 % 'local'
                keepMe(~r,~r) = 1;

        end
        % Remove missing data (NaNs encode no link):
        keepMe(isnan(allLinkData)) = 0;

        % Keep values assigned to each remaining link as this element of hubhubData
        linkDataSpecial = allLinkData(keepMe);
        allHubHub{1}{i,j} = linkDataSpecial;

        notSpecial = (~isnan(allLinkData) & ~keepMe);
        linkDataNotSpecial = allLinkData(notSpecial);

        % 2-sample t-test for special links greater than non-special links:
        if doRanksum
            if any(~isnan(linkDataSpecial)) && any(~isnan(linkDataNotSpecial))
                p = ranksum(linkDataSpecial,linkDataNotSpecial, 'tail', 'right');
            else
                p = NaN;
            end
        else
            [~,p] = ttest2(linkDataSpecial,linkDataNotSpecial,'Vartype','unequal', 'Tail','right');
        end
        tStats(i,j) = p;
    end
end

% ------------------------------------------------------------------------------
%% Plot as rich plots
% ------------------------------------------------------------------------------
myColors = GiveMeColors('richFeedInOutPeripheral'); % [BF_getcmap('spectral',4,1),BF_getcmap('set2',4,1)];
plotOnOne = true; % plot all on one figure
includeHist = true;
plotJustRich = true;
includeStd = false;
plotNulls = false; % plot each null trajectory in the figure
sigThresh = 0.05;

for j = 1:length(whatLinks)
    if plotJustRich && (j < 1)
        break
    end
    if plotDist
        JitteredParallelScatter(allHubHub{1}(:,j))
    else
        if includeHist
            if (~plotOnOne || (plotOnOne==1 && j==1))
                figure('color','w');
                sp=subplot(5,3,2:3);
                pos=get(sp,'Position');
                set(sp,'Position',[pos(1), pos(2)*0.94, pos(3), pos(4)]); % [left bottom width height]hold on
                % Degree distribution
                %N = arrayfun(@(x)sum(nodeData==x),krAll);
                %bar(krAll,N,'EdgeColor','k','FaceColor','k')
                histogram(nodeData,100,'EdgeColor','k','FaceColor','k');
                xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
                xticks([]); box off;
                ylabel('Frequency', 'FontSize', 18)
                get(gca, 'YTick');
                set(gca, 'FontSize', 16)
                % set(gca,'XTickLabel',{})
                sp=subplot(5,3,5:6);
                % % Set the Figure Size and Position (so that the labels fit)
                pos=get(sp,'Position');
                set(sp,'Position',[pos(1)*1, pos(2)*0.5, pos(3)*1, pos(4)*2.5]); % [left bottom width height]
                set(gca,'Ytick', [0 0.1 0.2 0.3 0.4 0.5], 'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5], 'FontSize', 18);
                hold on;

            end
        else
            if (~plotOnOne || (plotOnOne==1 && j==1))
                figure('color','w'); hold on
            end
        end
    end
    xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
    xlabel(extraParam{1})
    switch measure
        case 'mean'
            ylabel(sprintf('Mean %s ',analyzeWhat));
        case 'median'
            ylabel(sprintf('Median %s ',analyzeWhat));
    end


    % Plot flat line for rich
    if j==1
        switch measure
            case 'mean'
                plot([kr(1),krAll(end)],ones(2,1)*mean(allHubHub{1}{1,j}),':','color','k','LineWidth',3)
            case 'median'
                plot([kr(1),krAll(end)],ones(2,1)*median(allHubHub{1}{1,j}),':','color','k','LineWidth',3);
        end
    end
    switch measure
        case 'mean'
            realTrajectory = cellfun(@mean,allHubHub{1}(:,j));
        case 'median'
            realTrajectory = cellfun(@median,allHubHub{1}(:,j));
    end

    realStd = cellfun(@std,allHubHub{1}(:,j));

    % p-values from 2-sample t-test with unequal variances:
    pvalues = tStats(:,j);

    isSig = (pvalues < sigThresh); % significantly higher than null

    % mean (real data trajectory):
    lineStyle = '-'; markerStyle = 'o';

    if any(isSig)
        plot(kr(isSig),realTrajectory(isSig),markerStyle,'MarkerEdgeColor',myColors(j+1,:),...
            'MarkerFaceColor',brighten(myColors(j+1,:),+0.5),'LineWidth',1,'MarkerSize',9)
    end

    % mean trajectory:
    plot(kr,realTrajectory,lineStyle,'color',myColors(j+1,:),'LineWidth',2.5)

    % +/- std:
    if includeStd
        plot(kr,realTrajectory+realStd,lineStyle,'color',myColors(j+1,:),'LineWidth',1)
        plot(kr,realTrajectory-realStd,lineStyle,'color',myColors(j+1,:),'LineWidth',1)
    end

    xLimits = get(gca,'xlim'); yLimits = get(gca,'ylim');

    if ~plotJustRich
        text(xLimits(1)+0.1*diff(xLimits),yLimits(1)+0.9*diff(yLimits)-j/20,whatLinks{j},'color',myColors(j+1,:),'FontSize',18)
    end

    % Add proportion of nodes that are hubs:
    %     theYLim = ge    % Add proportion of nodes that are hubs:
    %     theYLim = get(gca,'ylim');
    %     plot(1:length(kr),theYLim(1)+(1-propisHub)*diff(theYLim),'--','color',myColors{4},'LineWidth',2)
    %
    %     % Add number of hub-hub links:
    %     plot(1:length(kr),theYLim(1) + (1-numHubHubLinks/max(numHubHubLinks))*diff(theYLim),'-.','color',brighten(myColors{4},-0.5),'LineWidth',2)
    %t(gca,'ylim');
    %     plot(1:length(kr),theYLim(1)+(1-propisHub)*diff(theYLim),'--','color',myColors{4},'LineWidth',2)
    %
    %     Add number of hub-hub links:
    %     plot(1:length(kr),theYLim(1) + (1-numHubHubLinks/max(numHubHubLinks))*diff(theYLim),'-.','color',brighten(myColors{4},-0.5),'LineWidth',2)

    % Add a meaningful title:
    divisionText = '';

    %     switch analyzeWhat
    %         case 'genecorr'
    %             %title(sprintf('%s correlations for genes across %s pairs of nodes of type %s (alpha = %.2g)', ...
    %             %networkType,whatNorm,whatType,sigThresh),'interpreter','none')
    %         otherwise
    %             title(sprintf('Mean %s analyzeWhat,sigThresh),...
    %                 'interpreter','none')
    %     end

end

ax = gca;
if plotJustRich
    % set(gcf,'Position',[1058,473,530,270])
    set(gcf,'Position',[1000,200,1400,700])
    ax.FontSize = 18;
else
    set(gcf,'Position',[1500,200,3000,2500])
end


% Add light gray rectangle
ylimNow = [ax.YLim(1),ax.YLim(2)];
if ax.YLim(1)<0
    h_rect = rectangle('Position',[44,ylimNow(1),20,ylimNow(2)+(ylimNow(1).*(-1))],'EdgeColor','none','FaceColor',ones(3,1)*0.90);
else
    h_rect = rectangle('Position',[44,ylimNow(1),20,ylimNow(2)-ylimNow(1)],'EdgeColor','none','FaceColor',ones(3,1)*0.90);
end
uistack(h_rect,'bottom');
set(gca,'ylim',ylimNow);

switch analyzeWhat
    case 'connection distance'
        %if pThreshold==0.05, ylim([0,1.2]); end
        ylabel('Average connection distance (mm)', 'FontSize', 18)
    case 'lineage distance'
        if pThreshold==0.05, ylim([14.5,18.5]); end
        ylabel('Average lineage distance', 'FontSize', 18)

        set(gca,'Ytick', [14.5 15 15.5 16 16.5 17 17.5 18 18.5], 'YTickLabel',[14.5 15 15.5 16 16.5 17 17.5 18 18.5]);
    case 'coexpression'
        switch coexpMeasure
            case 'coX'

                ylim([0.8 1]);
                set(gca,'Ytick', [0.8 0.85 0.9 0.95 1], 'YTickLabel',[0.8 0.85 0.9 0.95 1], 'FontSize', 18);
            case 'Pearson_noLR'
                ylim([0.0,0.5]);

                set(gca,'Ytick', [0 0.1 0.2 0.3 0.4 0.5], 'YTickLabel',[0 0.1 0.2 0.3 0.4 0.5], 'FontSize', 18);
                %set(gca,'Xtick', [0 10 20 30 40 50 60], 'YTickLabel',[0 10 20 30 40 50 60], 'FontSize', 14);
        end
        switch measure
            case 'mean'
                axisName = {'Mean gene', 'coexpression, r_\phi'};
            case 'median'
                axisName = {'Median gene', 'coexpression, r_\phi'};
        end
        ylabel(axisName, 'FontSize', 18)
    case 'birth time difference'
        ylim([-100 800]);
        set(gca,'Ytick', [0 200 400 600 800], 'YTickLabel',[0 200 400 600 800], 'FontSize',18);
end
%set(gca,'XTick',40:80)

xlabel('Degree, k','FontSize', 18); 

[link, deg] = propLinkDegree(C,whatAdj,D.whatConns);
% f = figure('color','w');
% subplot(,1,1)
% kRange = degsort;
% N = arrayfun(@(x)sum(deg==x),kRange);
% bar(kRange,N,'EdgeColor','k','FaceColor','k')
% xlim([min(deg),max(deg)]);
% ylabel('Frequency')
sp=subplot(5,3,5:6);
pos=get(sp,'Position');
set(sp,'Position',[pos(1), pos(2)*1.025, pos(3), pos(4)*0.7]); % [left bottom width height]
% Give colors:
%colors = BF_getcmap('spectral',4,1);
kRange = min(deg):max(deg);
for i = 1:length(kRange)
    %ind = (deg >= degsort(i));
    %categoriesHere = C.RegionM(ind==1, types);
    proportions = link'; %countcats(categoriesHere)/length(categoriesHere);
    for j = 1:size(proportions,2)
        if proportions(i,j) > 0
            rectangle('Position',[kRange(i)-1,sum(proportions(i,1:j-1)),1,proportions(i,j)], ...
                'FaceColor',myColors(j+1,:),'EdgeColor',myColors(j+1,:))
        end
    end
end
ylabel({'Proportion';'of links'},'FontSize', 18);
xlim([min(nodeData)-0.5,max(nodeData)+0.5]);
get(gca, 'YTick');
set(gca, 'FontSize', 16)
xticks([]);
ylim([0 1]);
yticks([0 1]);
hold on;
% add connected vs  unconnected plot
[dataCell, S,P] = coexpELCHUncon(C, G, true);
sp = subplot(5,3,7:7.5);
% % Set the Figure Size and Position (so that the labels fit)
pos = get(sp,'Position');
set(sp,'Position',[pos(1)*0.6, pos(2)*0.85, pos(3)*1.3, pos(4)*4]); % [left bottom width height]
hold on;
extraParams = struct('customSpot','.');
rgb_colorMatrix = GiveMeColors('ElChemUncon');
colors = num2cell(rgb_colorMatrix, 2);
extraParams.theColors = colors;

JitteredParallelScatter(dataCell,true,1,false,extraParams);
%L1 = {'Electrical';sprintf('(%d pairs)', length(dataCell{1}))};
set(gca,'Xtick', [1 2 3], 'XTickLabel',{sprintf('Electrical (%d pairs)', length(dataCell{1})), sprintf('Chemical (%d pairs)', length(dataCell{2})), sprintf('Unconnected (%d pairs)', length(dataCell{3}))}, 'FontSize', 18);
%xtickangle(30);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 18);
set(gca,'box','off');
ylabel('Gene coexpression, r_\phi','FontSize', 18);

end
