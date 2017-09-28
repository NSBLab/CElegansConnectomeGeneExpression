function probabilityVSdistance(C,numBins,doSymmetric,distinguishClasses,makeNewFigure)
% ------------------------------------------------------------------------------
% Function plots probability that a connection exists as a function of
% distance between neurons
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% numBins - choose number of bins for distance          default 20
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
if nargin < 2
    numBins = 50;
    fprintf(sprintf('%d bins by DEFAULT\n', numBins));
end
if nargin < 3
    doSymmetric = false;
end
if nargin < 4
    distinguishClasses = false;
end
if nargin < 5
    makeNewFigure = true;
end
%-------------------------------------------------------------------------------

Adj = GiveMeAdj(C,'zeroBinary');
dist = GiveMeDist(C);

if doSymmetric
    % Symmetrize to get probability two regions will be connected (at all)
    AdjSym = (Adj | Adj');
    % Get upper triangle as vectors:
    distVector = dist(triu(true(size(dist)),+1));
    connVector = AdjSym(triu(true(size(AdjSym)),+1));
else
    % Do not symmetrize to get probability that a connection exists for every
    % case in which a connection could exist
    % Get vectors as non-diagonal elements:
    notDiag = ~logical(eye(size(dist)));
    distVector = dist(notDiag);
    connVector = Adj(notDiag);
end

% Masks by type
if distinguishClasses
    maskAnatomy = GiveMeMask(C,'AnatomyPairs');
    [theColors,maskColorLabels] = GiveMeColors('directedAnatomy');

    if makeNewFigure
        f = figure('color','w'); hold on
        f.Position = [1000        1149         614         189];
    end
    maskNames = fieldnames(maskAnatomy);
    % hyphenate mask names for legend
    maskNamesHyphenated = cell(size(maskNames));
    for k = 1:length(maskNames)
        maskNamesHyphenated{k} = [maskNames{k}(1:4),'->',maskNames{k}(5:end)];
    end
    handles = cell(length(maskNames),1);

    groupings = cell(2,1);
    groupings{1} = {'HeadHead','TailTail','BodyBody','HeadTail','TailHead'};
    groupings{2} = {'HeadBody','BodyHead','BodyTail','TailBody'};
    ax = cell(2,1);
    for k = 1:2
        ax{k} = subplot(1,2,k);
    end
    whatGroup = zeros(length(maskNames),1);
    for k = 1:length(maskNames)
        maskAnatomy.(maskNames{k}) = logical(maskAnatomy.(maskNames{k})(notDiag));
        whatGroup(k) = find(cellfun(@(x)ismember(maskNames{k},x),groupings));
        dist_k = distVector(maskAnatomy.(maskNames{k}));
        conn_k = connVector(maskAnatomy.(maskNames{k}));
        subplot(1,2,whatGroup(k));
        hold on

        handles{k} = BF_PlotQuantiles(dist_k,conn_k,...
                                numBins,0,false,theColors(strcmp(maskColorLabels,maskNames{k}),:));
        % Fit exponential??
        % if ismember(maskNames{i},{'HeadHead','TailTail'})
        xData = dist_k;
        yData = conn_k;
        if ismember(maskNames{k},{'HeadTail','TailHead'})
            s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,1,0]);
            f = fittype('A*exp(-n*(x - C)) + B','options',s);
            hasOffset = true;
            beep
        else
            s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,0]);
            f = fittype('A*exp(-n*x) + B','options',s);
            hasOffset = false;
        end
        try
            [c, Stats] = fit(xData,yData,f);
        catch
            error('Fit to exp failed')
        end
        if hasOffset
            f_handle = @(x) c.A.*exp(-c.n*(x-c.C)) + c.B;
        else
            f_handle = @(x) c.A.*exp(-c.n*x) + c.B;
        end
        % ---
        xRange = linspace(min(xData),max(xData),100);
        plot(xRange,f_handle(xRange),':','color',theColors(strcmp(maskColorLabels,maskNames{k}),:),'LineWidth',2)
        % end
    end
    for k = 1:2
        legend(ax{k},[handles{whatGroup==k}],maskNamesHyphenated(whatGroup==k),'Location','North');
        ax{k}.YLim = [0,0.2];
        ax{k}.XLim = [0,1.2];
        if k==1
            LabelCurrentAxes('A',ax{k},18,'topRight');
        else
            LabelCurrentAxes('B',ax{k},18,'topRight');
        end
    end

    % Put axes together:
    ylabel(ax{1},'Connection probability');
    xlabel(ax{1},'Separation distance (mm)');
    ylabel(ax{2},'Connection probability');
    xlabel(ax{2},'Separation distance (mm)');
    % linkaxes([ax{1},ax{2}],'x')
else
    BF_PlotQuantiles(distVector,connVector,numBins,0,makeNewFigure);
    xlabel('Separation distance (mm)');
    ylabel('Connection probability');
end

end
