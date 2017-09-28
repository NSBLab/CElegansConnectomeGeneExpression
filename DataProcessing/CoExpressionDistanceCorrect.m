function coExp_corr = CoExpressionDistanceCorrect(C,G)
%-------------------------------------------------------------------------------
% Corrects for distance dependence in coexpression
% Takes in:
% G -> pairwise coexpression matrix, and
% C -> pairwise Euclidean distance matrix
%-------------------------------------------------------------------------------

distMat = GiveMeDist(C);
coExpMat = GiveMeCoexpression(G,'',true);
fprintf(1,'Using Pearson for coexpression\n');
coExp_corr = coExpMat;

% Only look at upper triangle
fprintf(1,'Looking at all pairs, disregarding axonal connectivity\n');
maskBasic = triu(true(size(distMat)),+1); % upper-triangle

maskAnatomy = GiveMeMask(C,'AnatomyPairs');

anatomyCorrect = {'HeadHead','TailTail'};

%-------------------------------------------------------------------------------
% Correct head-head and tail-tail pairs
%-------------------------------------------------------------------------------
for i = 1:length(anatomyCorrect)
    combinedMask = maskBasic & maskAnatomy.(anatomyCorrect{i});

    xData = distMat(combinedMask);
    yData = coExpMat(combinedMask);
    isGood  = ~isnan(yData);
    xData = xData(isGood);
    yData = yData(isGood);

    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,0.5,0]);
    f = fittype('A*exp(-n*x) + B','options',s);
    try
        [c, Stats] = fit(xData,yData,f);
        f_handle = @(x) c.A.*exp(-c.n*x) + c.B;
        % ---
        % xRange = linspace(min(xData),max(xData),100);
        % plot(xRange,f_handle(xRange),'color',anatomyColors(j,:),'LineWidth',2.5)
    catch
        error('Fit to exp failed')
    end

    % Correct to residuals:
    coExp_corr(combinedMask) = coExpMat(combinedMask) - f_handle(distMat(combinedMask));

    % Plot:
    f = figure('color','w');
    f.Position = [162   859   519   214];
    ax1 = subplot(1,2,1); hold on
    plot(distMat(combinedMask),coExpMat(combinedMask),'.','color',ones(1,3)*0.5)
    xRange = linspace(min(xData),max(xData),100);
    xlabel('Separation distance (mm)')
    ylabel('Correlated gene expression, r_\phi')
    title(anatomyCorrect{i})
    plot(xRange,f_handle(xRange),'k','LineWidth',2.5)
    LabelCurrentAxes('A',gca,18,'topRight');
    ax2 = subplot(1,2,2);
    plot(distMat(combinedMask),coExp_corr(combinedMask),'.','color',ones(1,3)*0.5)
    xlabel('Separation distance (mm)')
    ylabel('Correlated gene expression residuals')
    box('off')
    LabelCurrentAxes('B',gca,18,'topRight');
    axYLims = [min([ax1.YLim(1),ax2.YLim(1)]),max([ax1.YLim(2),ax2.YLim(2)])];
    ax1.YLim = axYLims;
    ax2.YLim = axYLims;
end

end
