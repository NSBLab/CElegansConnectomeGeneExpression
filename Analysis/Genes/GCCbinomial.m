function [T,geneScores] = GCCbinomial(C,G,chooseType,countON,DOsave)
%-------------------------------------------------------------------------------
% Aurina Arnatkeviciute 23-01-2017
% Aurina Arnatkeviciute update 09-02-2017
% Idea is to provide a score for each gene according to the number of
% matches in gene expression "both on" or "both on/both off" falling on
% specific link type. Will focus on connected links and rich&feeder links.
% countON - 1 - count only "both on" mathes. Any other number would give both "on&OFF" matches.
% khub - define degree for hub threshold
% chooseType - 'Connected', 'RichFeeder'
%-------------------------------------------------------------------------------
% load connectivity matrix
if nargin < 3
    chooseType = 'Connected';
    fprintf('Connected by DEFAULT\n')
    % can choose RichFeeder
end

if nargin < 4
    countON = 1;
    fprintf('Count ON matches only by DEFAULT\n')
end

if nargin < 5
    DOsave = false;
    fprintf('NOT saving output by DEFAULT\n')
end

thresholdBadMatches = 10; % if this many matches (or fewer), exclude.
removeBilateral = true; % exclude possibility of bilateral matches

%-------------------------------------------------------------------------------
Adj = GiveMeAdj(C,'zeroBinary');
D = GiveMeDefault();
kHub = D.kHub;
masks = GiveMeMask(C,'hubType','',Adj);
connMask = (Adj|Adj'); % symmetrize
connMask(tril(true(279))) = 0; % ignore lower diagonal


numNeurons = C.numNeurons;
numGenes = size(G.GeneExpData.Direct,2);

% Calculate the number of connected links in the connectome
%numConnected = sum(connMask(:)); % number of connected neuron pairs
%numPossiblePairs = numNeurons*(numNeurons-1)/2; % number of possible neuron pairs

%-------------------------------------------------------------------------------
% Settings
%-------------------------------------------------------------------------------
% Define: theMask -- the set of pairs that we're looking for matches on
% theBase -- the base set of pairs where matches can be
switch chooseType
    case 'Connected'
        theMask = connMask; % on connections
        theBase = triu(true(279),+1); % all pairs
    case 'Rich'
        theMask = masks.rich; % on hub-hub
        theBase = connMask; % on connected pair
    case 'RichFeeder'
        theMask = masks.richfeed;
        theBase = connMask;
    case 'RichFeedIN'
        theMask = masks.richfeedin;
        theBase = connMask;
    case 'RichFeedOUT'
        theMask = masks.richfeedout;
        theBase = connMask;
    case 'FeedIN'
        theMask = masks.feedin;
        theBase = connMask;
    case 'FeedOUT'
        theMask = masks.feedout;
        theBase = connMask;
    case 'Chemical'
        theMask = masks.chemicalOnly;
        theBase = masks.CHandELonly;
    case 'Electrical'
        theMask = masks.electricalOnly;
        theBase = masks.CHandELonly;
end

%-------------------------------------------------------------------------------
% Make sure no bilateral components can exist
if removeBilateral
    fprintf(1,'~~~Bilateral pairs are not considered in the analysis!!!\n');
    maskIsLR = GiveMeLRMask(C);
    theMask(maskIsLR) = 0;
    theBase(maskIsLR) = 0;
end

matchProb = sum(theMask(:))/sum(theBase(:)); % probability of being on the mask, given on the base

%===============================================================================
% Loop through all genes and calculate the score
%===============================================================================
matchesMask = zeros(numGenes,1);
matchesBase = zeros(numGenes,1);
geneScore = zeros(numGenes,1);
for i = 1:numGenes
    % get numbers of matching on and off instances for a real gene
    gene = G.GeneExpData.Direct(:,i);
    match1 = gene*gene';
    match0 = (1-gene)*(1-gene)';

    % Define the matching matrix, M: 1 for gene expression matches
    switch countON
    case 0
        M = triu(match1|match0,1);
        fileN = 'ONOFF';
    case 1
        M = triu(match1,1); %|match0;
        fileN = 'ON';
    case 2
        M = triu(match0,1);
        fileN = 'OFF';
    end
    M = double(M);
    % M(logical(eye(size(M)))) = NaN; % exclude diagonal

    %-------------------------------------------------------------------------------
    % Compute the number of matches on the mask
    M_base = M.*theBase;
    M_mask = M.*theMask;
    matchesBase(i) = sum(M_base(:));
    matchesMask(i) = sum(M_mask(:));

    % Compute the probability of obtaining more matches than observed:
    if matchesBase(i) <= thresholdBadMatches % if there are too few matches
        % if there one or less matches falling on mask (here - on connected links) we can't say anything about matches on rich/feeder links.
        geneScore(i) = NaN;  % excluding these genes makes a lot of difference in enrichment as comparison is performed within a subset and we get rid of around half genes. Enrichment fo rrich and feeder links loose the signifficance.
    else
         if matchesMask(i)==0
             geneScore(i) = 1;
         else
            geneScore(i) = 1 - binocdf(matchesMask(i)-1,matchesBase(i),matchProb);
        end
    end

    % if we want to take -log(score), then will need to assign some
    % value to 0 (as they turn as Inf when loged).
    % ermineJ gives zeros --> 10^-15 value.
end

%-------------------------------------------------------------------------------
% do FDR correction

% Get gene names
geneName = G.geneAcronyms.Direct;

% Exclude NaN genes
isBadGene = isnan(geneScore);
fprintf(1,'%u/%u genes had too few matches...\n',sum(isBadGene),length(isBadGene));
geneName(isBadGene) = [];
matchesMask(isBadGene) = [];
matchesBase(isBadGene) = [];
geneScore(isBadGene) = [];
% do FDR correction
geneScores_corr = mafdr(geneScore,'BHFDR',true);

% Output to table:
T = table(geneName,geneScore,matchesMask,matchesBase,geneScores_corr);
geneScores = table(geneName,geneScores_corr);

if DOsave
    if removeBilateral
        extraText = 'NOLR';
    else
        extraText = '';
    end
    fileNameOut = fullfile('Data','ermineJdata',sprintf('Corrected%sscores_%s_k%d%sN%dATLEASTn.txt',chooseType,fileN,kHub,extraText,thresholdBadMatches));
    writetable(geneScores,fileNameOut,'Delimiter','\t')
end

end
