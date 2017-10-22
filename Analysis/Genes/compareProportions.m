function [BinCorrelation,averageCorrelation] = compareProportions(binaryCorrelation,makeFigure)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-11-17
% This script compares how different correlation measures depend on the
% proportion of on and off instances in binary variables
% ------------------------------------------------------------------------------

if nargin < 2
    makeFigure = true;
end

minProp = 1; % min number of ones in a vector
maxProp = 148; % max number of ones in a vector - proportion for our data is 0.144
numGenes = 948; % length(G.GeneStruct);
naverage = 1000; % how many runs to average over
numValues = 148; % number of vectors to calculate coexpression on.

% create a random vector with a set proportion of 1;
numCases = round(linspace(minProp,maxProp,numValues));
M1 = zeros(length(numCases),naverage,numGenes);
for j=1:length(numCases)
    for l=1:naverage
        m = numCases(j);
        randVect = zeros(numGenes,1);
        randVect(randperm(numGenes,m)) = 1;
        %allrandVect (l,:) = randVect;
        M1(j,l,:) = randVect;
    end
end

%M1 = squeeze(sampleVectors(:,1,:));
%M2 = squeeze(sampleVectors(:,2,:));

% calculate correlations between those vectors
numVectors = length(numCases);
BinCorrelation = zeros(numVectors,numVectors,naverage);
switch binaryCorrelation
case 'coX'
    for k=1:naverage
        fprintf(1,'Iteration %u/%u\n',k,naverage);
        M2 = squeeze(M1(:,k,:));
        BinCorrelation(:,:,k) = CoExpressionPositiveMatch(M2,true,true,false);
    end
otherwise
    for k=1:naverage
        fprintf(1,'Iteration %u/%u\n',k,naverage);
        for i = 1:numVectors
            for j = 1:numVectors
                BinCorrelation(i,j,k) = BinaryCorr(squeeze(M1(i,k,:)),squeeze(M1(j,k,:)),binaryCorrelation);
            end
        end
    end
end
%BinCorrelation(1:size(BinCorrelation,1)+1:end) = nan;
averageCorrelation = mean(BinCorrelation,3);
averageCorrelation(1:size(averageCorrelation,1)+1:end) = nan;

if makeFigure
    figure('color','w');
end
ax = gca;
imagesc(averageCorrelation);
title(sprintf('Correlation matrix %s', binaryCorrelation));
colorbar; %caxis([-1 1]);
colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]); %
xlabel(sprintf('Number of ones in a vector 1 (/%u)',numGenes));
ylabel(sprintf('Number of ones in a vector 2 (/%u)',numGenes));
xtick = numCases;
ax.XTick = 1:length(numCases);
ax.XTickLabel = xtick;
ax.YTick = ax.XTick;
ax.YTickLabel = ax.XTickLabel;
% ([1,10:10:150]); xtick = round(100*xtick/numGenes)/100;
% set(gca,'XTick',[1,10:10:150],'XTickLabel',xtick);
% set(gca,'YTick',[1,10:10:150],'YTickLabel',xtick);

%-------------------------------------------------------------------------------
    function rho = BinaryCorr(x,y,CorrMethod)
        % 1. Calculate contingency table
        a = sum(x==1&y==1);
        b = sum(x==1&y==0);
        c = sum(x==0&y==1);
        d = sum(x==0&y==0);
        n = a+b+c+d;
        p1 = a+b;
        q1 = c+d;
        p2 = a+c;
        q2 = b+d;
        switch CorrMethod
            case 'Loevinger'
                rho = (a*d-c*b)/(min(p1*q2, p2*q1));
            case 'Pearson'
                % 2. Calculate correlation coefficient
                rho = (a*d-c*b)/sqrt((a+b)*(c+d)*(a+c)*(b+d));
            case 'SimpleMatching'
                rho = (a+d)/(a+b+c+d);
            case 'Yule'
                rho = ((a*d)-(b*c))/((a*d)+(b*c));
            case 'YuleY'
                rho = (sqrt(a*d)-sqrt(b*c))/(sqrt(a*d)+sqrt(b*c));
            case 'Jaccard'
                rho = a/(a+b+c);
            case 'Dice'
                rho = (2*a)/((2*a)+b+c);
            case 'Prop'
                rho = a/min((a+b),(a+c));
            case 'Chi'
                rho = sqrt((n*(a*d - c*b)^2)/((a+b)*(c+d)*(a+c)*(b+d)));
        end
    end
end
