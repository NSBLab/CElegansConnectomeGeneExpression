function pValues = PlotRichClubLinked(C,networkType,whatType,sym,pAdjInclude,includeHist)
% Plots a pre-computed set of conventional rich club coefficients
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-10-08 -- turned to a function
% Modified by Taqi Ali 17-7-15
%
% sym{1,0} = whether the adjacency matrix used to find data is symmetric or not
% (applies to Towlsons undirected form of the "all" network which includes
% Ch and EJ junctions)
% Other variables: See possible alternatives
% ------------------------------------------------------------------------------

% Analyze either a binary or weighted (directed) network
if nargin < 2
    networkType = 'bd'; %'bu','bd','wd','wu'
end
if nargin < 3
    whatType = 'all'; %'ch','ej','all'
end
if nargin < 4
    sym = 0;
end
if nargin < 5
    pAdjInclude = 0.05;
end
if nargin < 6
    includeHist = 1;
end

if sym == 1;
    issym = 'sym';
else
    issym = 'non-sym';
end


plotStyle = 'single'; % 'single' (single figure), 'multiple' (subplots in figure).


pThreshold = 0.05; % FOR SIGNIFICANT PHI
numComparisons = 1; % number of unique degrees (59)

% ------------------------------------------------------------------------------
% Load data:
% ------------------------------------------------------------------------------
fileName = sprintf('RichClub_%s_%s_%s_%.0f.mat',networkType,whatType,issym,pAdjInclude*100);
load(fileName);
fprintf(1,'Loaded data from %s with %u sets of randomized rich club data computed\n',...
                                    networkType,size(PhiRand,1));
C = load('CElegans_Connectivity_Data.mat');

% ------------------------------------------------------------------------------
% Compute Phi:
% ------------------------------------------------------------------------------

% Get individual phi_norm wrt different randomizations:
PhiNorm = zeros(numRepeats,kmax);
for i = 1:numRepeats
    PhiNorm(i,:) = PhiTrue./PhiRand(i,:);
end

% Ben Fulcher, 2014-11-13 (error found in the process of giving feedback to Simon Baker)
% Recompute PhiNormMean as the ratio of the true Phi, to the mean across randomizations:
PhiNormMean = zeros(size(PhiTrue));
for i = 1:length(PhiTrue)
    PhiNormMean(i) = PhiTrue(i)/mean(PhiRand(:,i));
end

% Compute p-values
pValues = zeros(kmax,1);
for i = 1:kmax
    pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
end

% Significant phi
isSig = (pValues <= pThreshold/numComparisons);

% ------------------------------------------------------------------------------
% Plotting:
myColors = BF_getcmap('spectral',4,1);
f = figure('color','w'); hold on;

switch plotStyle
case 'single'
    if includeHist
        kr = min(k_tot):max(k_tot);
        linkedAdj = logical(GiveMeAdj(C,'rawBinary',whatType));
        switch networkType
            case {'bu', 'wu'}
                nodeData = sum(linkedAdj);
            case {'bd', 'wd'}
                nodeData = sum(linkedAdj) + sum(linkedAdj');
        end
        subplot(6,1,1:2); hold on;
        % Degree distribution:
        N = arrayfun(@(x)sum(nodeData==x),kr);
        bar(kr,N,'EdgeColor','k','FaceColor','k'); %myColors{2})
        xlim([min(k_tot)-0.8,max(k_tot)+0.8]);
        ylabel('Frequency');

        title(sprintf(['Rich club curves for %s networks, of type %s, using %u ' ...
                        'iterations, for %u different nulls and method %s'], ...
                    networkType,whatType,numIter,numRepeats,nullMethod),'interpreter','none');
        ax = gca;
        ax.XTickLabel = {};
        ax.FontSize = 11;

        subplot(6,1,3:6); hold on;
    end

    % Plot a single figure: the phi_norm curves
    plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor',myColors{1},...
                        'MarkerFaceColor',brighten(myColors{1},+0.5),'MarkerSize',6,'LineWidth',1)
    plot(PhiNormMean,'-','color',myColors{1},'LineWidth',2);
    plot([1,kmax],ones(2,1),':k');
    ylabel('\Phi_{norm}');
    xlabel('Degree, k');

    xlim([min(k_tot)-0.8,max(k_tot)+0.8]);

    % Reposition the figure:
    if includeHist
        % f.Position = [458,557,386,276]; % bigger
        f.Position = [458,557,386,207]; % smaller
    else
        f.Position = [458,557,386,207];
        % f.Position = [1058,473,530,270];
    end

    % Add light gray rectangle:
    if pAdjInclude==0.05
        ax = gca;
        h_rect = rectangle('Position',[42,ax.YLim(1),18,ax.YLim(2)-ax.YLim(1)],'EdgeColor','none','FaceColor',ones(3,1)*0.90);
        uistack(h_rect,'bottom');
        ax.FontSize = 11;
    end

case 'multiple'
    % Plot four sub-figures
    subplot(2,2,1); hold on; % the phi_norm curves
    plot(PhiNormMean,'-','color',myColors{1},'LineWidth',2);
    plot(find(isSig),PhiNormMean(isSig),'o','color',myColors{1},'LineWidth',2);
    plot([1,kmax],ones(2,1),'-.k');
    ylabel('\phi_{norm}');
    xlabel('Degree');
    xlim([min(k_tot),max(k_tot)]);
    title(sprintf(['Rich club curves for %s networks, of type %s, using %u ' ...
                        'iterations, for %u different nulls and method %s'], ...
                    networkType,whatType,numIter,numRepeats,nullMethod),'interpreter','none');

    subplot(2,2,2); hold on; % the phi curves
    h1 = plot(mean(PhiRand),'o-','color',ones(3,1)*0.4,'LineWidth',2);
    h2 = plot(PhiTrue,'-','color',myColors{1},'LineWidth',2);
    h3 = plot(find(isSig),PhiTrue(isSig),'o','color',myColors{1},'LineWidth',2);
    h4 = plot(PhiNormMean,'-','color',myColors{2},'LineWidth',2);
    xlim([min(k_tot),max(k_tot)]);
    ylabel('\phi');
    xlabel('Degree');
    title(['Phi curves for ' networkType ' network of type ' whatType],'interpreter','none');
    legend([h1 h2 h4],{'PhiRand','PhiTrue','PhiNor'},'Location','southeast');
    subplot(2,2,3); hold on; % the phi_norm curves
    plot(PhiNorm',':k');
    plot(PhiNormMean,'-','color',myColors{1},'LineWidth',2);
    plot(find(isSig),PhiNormMean(isSig),'o','color',myColors{1},'LineWidth',2);
    plot([1,kmax],ones(2,1),'-.k');
    ylabel('\phi_{norm}');
    xlabel('Degree');

    subplot(2,2,4); hold on; % the phi curves
    ylabel('\phi');
    xlabel('Degree');
    plot(PhiRand',':k');
    plot(mean(PhiRand),'o-','color',ones(3,1)*0.4,'LineWidth',2);
    plot(PhiTrue,'o-','color',myColors{1},'LineWidth',2);
end
    figure
    hold on
    h1 = plot(mean(PhiRand),'o-','color',ones(3,1)*0.4,'LineWidth',2);
    h2 = plot(PhiTrue,'-','color',myColors{1},'LineWidth',2);
    h3 = plot(find(isSig),PhiTrue(isSig),'o','color',myColors{1},'LineWidth',2);
    h4 = plot(PhiNormMean,'-','color',myColors{2},'LineWidth',2);
    xlim([min(k_tot),max(k_tot)])
    ylabel('\phi')
    xlabel('Degree')
    title(['Phi curves for ' networkType ' network of type ' whatType],'interpreter','none')
    legend([h1 h2 h4],{'PhiRand','PhiTrue','PhiNor'},'Location','southeast');
end
