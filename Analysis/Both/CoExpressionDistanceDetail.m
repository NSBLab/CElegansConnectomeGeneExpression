function CoExpressionDistanceDetail(C,G)
% Quick script to look into detail of coexpression relationship (within head-head)
% (or tail-tail)
% ------------------------------------------------------------------------------
doConnected = false;
makeNewFigure = true;
%-------------------------------------------------------------------------------
% Get pairwise distances, and coexpression values
distMat = GiveMeDist(C);
coExpMat = GiveMeCoexpression(G,[],true); 

% Only keep connections
if doConnected
    % Binary chemical adjacency matrix:
    Adj = GiveMeAdj(C,'zeroBinary');
    fprintf(1,'Only looking at where connections exist\n');
    mask = Adj & ~logical(eye(size(Adj))); % off-diagonal connections
else
    % Only look at upper triangle
    fprintf(1,'Looking at all pairs, disregarding axonal connectivity\n');
    mask = triu(true(size(distMat)),+1); % upper-triangle
end

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
lineWidth = 1.5;
numBins = 8;
addStd = true;
alsoScatter = false;

maskAnatomy = GiveMeMask(C,'AnatomyPairs');
[anatomyColors,maskLabels] = GiveMeColors('directedAnatomy');
if doConnected
    maskType = GiveMeMask(C,'TypePairs');
else
    maskType = GiveMeMask(C,'TypePairsSymmetric');
end
typeNames = fieldnames(maskType);
numTypes = length(typeNames);

anatomyDetail = {'HeadHead','TailTail'};

for i = 1:2
    f = figure('color','w');
    f.Position = [234   660   724   434];
    combinedMask = mask & maskAnatomy.(anatomyDetail{i});
    % Scatter (by type):
    handles = cell(numTypes,1);
    for j = 1:numTypes
        if doConnected
            ax = subplot(3,3,j); hold on
        else
            ax = subplot(2,3,j); hold on
            if j > 3
                xlabel('Separation distance (mm)')
            end
            if ismember(j,[1,3,5])
                ylabel('Gene coexpression, r_\phi')
            end
        end
        maskTot = combinedMask & maskType.(typeNames{j});
        distVec_i = distMat(maskTot);
        coExpVec_i = coExpMat(maskTot);
        isGood = ~isnan(coExpVec_i);
        if ~all(isGood)
            distVec_i = distVec_i(isGood);
            coExpVec_i = coExpVec_i(isGood);
            fprintf(1,'%u NaNs removed from coexpression for %s-%s\n',sum(~isGood),anatomyDetail{i},typeNames{j});
        end
        numConns = sum(maskTot(:));
        if numConns > 100
            theMarker = '.';
        else
            theMarker = 'x';
        end
        plot(distVec_i,coExpVec_i,theMarker,'color',brighten(anatomyColors(j,:),0.5));
        BF_PlotQuantiles(distMat(maskTot),coExpMat(maskTot),...
                                numBins,0,false,anatomyColors(j,:));
        % Fit exponential??
        xData = distVec_i;
        yData = coExpVec_i;
        s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,0]);
        f = fittype('A*exp(-n*x) + B','options',s);
        try
            [c, Stats] = fit(xData,yData,f);
            f_handle = @(x) c.A.*exp(-c.n*x) + c.B;
            % ---
            xRange = linspace(min(xData),max(xData),100);
            plot(xRange,f_handle(xRange),'color',anatomyColors(j,:),'LineWidth',2.5)
        catch
            warning('Fit to exp failed')
        end
        title(sprintf('%s-%s (%u)',anatomyDetail{i},typeNames{j},numConns))
        % Label axes A--F:
        switch j
        case 1, t = 'A';
        case 2, t = 'B';
        case 3, t = 'C';
        case 4, t = 'D';
        case 5, t = 'E';
        case 6, t = 'F';
        end
        LabelCurrentAxes(t,ax,16,'topLeft');
    end
    % handles{i} = BF_PlotQuantiles(distVec_i,coExpVec_i,numBins,alsoScatter,false,anatomyColors(whatColorIndex,:),addStd);
end

end
