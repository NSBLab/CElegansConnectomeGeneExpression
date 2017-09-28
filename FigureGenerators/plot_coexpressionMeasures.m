%===============================================================================
% SUPPLEMENTARY TEXT:
%===============================================================================
% Plot how different coexpression measures depend on the proportion of ones
% and zeros in a vector
%-------------------------------------------------------------------------------

f = figure('color','w');

measuresToCompare = {'Pearson','Yule','Jaccard','coX'};
for i = 1:4
    subplot(2,2,i)
    compareProportions(measuresToCompare{i},false);
    axis('square')
end

%===============================================================================
% Load in data for possitive correlation measure (too long to compute on (PC)
%===============================================================================

% load('propPositive.mat');
% minProp = 1; % min number of ones in a vector
% maxProp = 150; % max number of ones in a vector - proportion for our data is 0.144
% numGenes = 948; % length(G.GeneStruct);
% naverage = 10; % how many runs to average over
% numValues = 150; % number of vectors to calculate coexpression on.
%
% numCases = round(linspace(minProp,maxProp,numValues));
% averageCorrelation = squeeze(mean(coexp,1));
% averageCorrelation(1:size(averageCorrelation,1)+1:end) = nan;
%
% figure('color','w');
% imagesc(averageCorrelation);
% title('Correlation matrix positive matches'); colorbar; %caxis([-1 1]);
% colormap([flipud(BF_getcmap('blues',9));[1 1 1];BF_getcmap('reds',9)]);
% xlabel('Proportion of ones in a vector 1'); ylabel('Proportion of ones in a vector 2')
% xtick = numCases([1 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]);
% xtick = round(100*xtick/numGenes)/100;
% set(gca,'Xtick', [1,10:10:150], 'XTickLabel',xtick);
% set(gca,'Ytick', [1,10:10:150], 'YTickLabel',xtick);
% set(gcf,'color','w');
