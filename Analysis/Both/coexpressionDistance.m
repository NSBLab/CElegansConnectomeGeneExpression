function coexpressionDistance(C,G,doConnected,coexpMeasure,makeNewFigure)
% ------------------------------------------------------------------------------
% function plots coexpression as a dunction of distance for C elegans data.
% head-head, tail-tail, and all other links are repressented in different colors.
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure), G (gene structure) should be loaded;
% doConnected: true - plot only coexpression for connected links, false -
% plot coexpressio for all possible neuron pairs.
% coexpMeasure - choose coexpression measure as defined in G.Corr
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-12-20
% ------------------------------------------------------------------------------

D = GiveMeDefault();
if nargin < 3
    doConnected = false;
end
if nargin >= 4 && ~isempty(coexpMeasure)
    D.coexpMeasure = coexpMeasure;
end
if nargin < 5
    makeNewFigure = true;
end

%-------------------------------------------------------------------------------
% Get pairwise distances, and coexpression values

distMat = GiveMeDist(C);
coExpMat = GiveMeCoexpression(G,'',true);

% Only keep connections
if doConnected
    % Make binary chemical adjacency matrix:
    Adj = GiveMeAdj(C,'zeroBinary');
    % Adj = double(Adj|Adj'); % symmetrize
    fprintf(1,'Only looking at where connections exist\n');
    mask = Adj & ~logical(eye(size(Adj))); % off-diagonal connections
else
    % Only look at upper triangle
    fprintf(1,'Looking at all pairs, disregarding axonal connectivity\n');
    mask = triu(true(size(distMat)),+1); % upper-triangle
end

% Mask to filter to relevant parts -> vectors
distVec = distMat(mask);
coExpVec = coExpMat(mask);

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
if makeNewFigure
    hFig = figure('color','w');
    hFig.Position = [94   687   818   211]; %[376   625   248   461];
end
hold on;
markerSize = 10;
lineWidth = 1.5;
numBins = 8;
addStd = true;
alsoScatter = false;

% sc{1} = scatter(distVec,coExpVec,markerSize,'MarkerEdgeColor',[0 .5 .5],...
%     'MarkerFaceColor',[0 .7 .7],'LineWidth',lineWidth);

if doConnected
    maskAnatomy = GiveMeMask(C,'AnatomyPairs');
else
    maskAnatomy = GiveMeMask(C,'AnatomyPairsSymmetric');
end

[anatomyColors,maskLabels] = GiveMeColors('directedAnatomy');
groupings = cell(3,1);
groupings{1} = {'HeadHead'};
groupings{2} = {'TailTail'};
groupings{3} = {'BodyBody','HeadBody','BodyHead','BodyTail','TailBody','HeadTail','TailHead'}; % inter-body
doScatter = {'HeadHead','TailTail'};
numGroupings = length(groupings);

ax = cell(numGroupings,1);
for k = 1:numGroupings
    ax{k} = subplot(1,numGroupings,k); hold on;
    if k < 3, ax{k}.YLim = [-0.2,1]; end
    if k==3, ax{k}.XLim = [0,1.2]; end
end
maskNames = fieldnames(maskAnatomy);
numMasks = length(maskNames);
whatGroup = zeros(numMasks,1);
for i = 1:numMasks
    combinedMask = mask & maskAnatomy.(maskNames{i});
    distVec_i = distMat(combinedMask);
    coExpVec_i = coExpMat(combinedMask);
    isGood = ~isnan(coExpVec_i);
    if ~all(isGood)
        distVec_i = distVec_i(isGood);
        coExpVec_i = coExpVec_i(isGood);
        fprintf(1,'%u NaNs removed from coexpression for %s\n',sum(~isGood),maskNames{i});
    end
    whatColorIndex = strcmp(maskLabels,maskNames{i});
    whatGroup(i) = find(cellfun(@(x)ismember(maskNames{i},x),groupings));

    subplot(1,numGroupings,whatGroup(i));
    % sc{i} = scatter(distVec_i,coExpVec_i,markerSize,'MarkerEdgeColor',anatomyColors(whatColorIndex,:),...
    %     'MarkerFaceColor',anatomyColors(whatColorIndex,:),...
    %     'LineWidth',lineWidth);
    numThresholds = min([round(length(distVec_i)/10),8]);
    alsoScatter = ismember(maskNames{i},doScatter);
    handles{i} = BF_PlotQuantiles(distVec_i,coExpVec_i,numThresholds,alsoScatter,false,anatomyColors(whatColorIndex,:),addStd);
    % Fit exponential??
    if ismember(maskNames{i},{'HeadHead','TailTail'})
        xData = distVec_i;
        yData = coExpVec_i;
        s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,0]);
        f = fittype('A*exp(-n*x) + B','options',s);
        try
            [c, Stats] = fit(xData,yData,f);
        catch
            error('Fit to exp failed')
        end
        f_handle = @(x) c.A.*exp(-c.n*x) + c.B;
        % ---
        xRange = linspace(min(xData),max(xData),100);
        plot(xRange,f_handle(xRange),'color',anatomyColors(whatColorIndex,:),'LineWidth',3)
    end
end

%-------------------------------------------------------------------------------
% Axes tweaks/labels
%-------------------------------------------------------------------------------
yLimFullMin = min(cellfun(@(x)x.YLim(1),ax));
yLimFullMax = max(cellfun(@(x)x.YLim(2),ax));
% hyphenate mask names for legend
maskNamesHyphenated = cell(size(maskNames));
for k = 1:length(maskNames)
    maskNamesHyphenated{k} = [maskNames{k}(1:4),'-',maskNames{k}(5:end)];
end
for k = 1:numGroupings
    % pause(1)
    ylabel(ax{k},'Correlated gene expression, r_\phi');
    xlabel(ax{k},'Separation distance (mm)');
    ax{k}.YLim = [yLimFullMin,yLimFullMax];
    legend(ax{k},[handles{whatGroup==k}],maskNamesHyphenated(whatGroup==k));
    switch k
    case 1
        LabelCurrentAxes('A',ax{k},16,'topLeft');
    case 2
        LabelCurrentAxes('B',ax{k},16,'topLeft');
    case 3
        LabelCurrentAxes('C',ax{k},16,'topLeft');
    end
end

end
