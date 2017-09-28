function [probMatches,numOnes,numMatches] = preComputePositiveMatch(G,doExact,whatMax,justUnique)
%-------------------------------------------------------------------------------
% Idea is that you take in gene expression data and precompute all calculations
% needed for the coexpression positive match measure
%-------------------------------------------------------------------------------

if nargin < 2
    doExact = true;
end
if nargin < 3
    whatMax = 10; % only up to a max of 10 matches
end
if nargin < 4
    justUnique = false;
end
%-------------------------------------------------------------------------------

gData = [G.GeneExpData.Direct];
numGenes = size(gData,2);
if justUnique
    numOnes = unique(sum(gData));
else
    numOnes = 1:max(sum(gData));
end
numUnique = length(numOnes);

% For every pairwise combination of the number of ones, we want to compute
% the probability of matches up to the maximum

if strcmp(whatMax,'max')
    numMatches = 0:max(numOnes);
else
    numMatches = 0:whatMax;
end

probMatches = nan(numUnique,numUnique,length(numMatches));
for i = 1:numUnique
    fprintf(1,'%u/%u: ',i,numUnique);
    for j = i:numUnique % i<=j
        fprintf(1,'%u,',j);
        for k = 1:min([max(numMatches),i])+1
            probMatches(i,j,k) = probMMatches(numGenes,numOnes(j),numOnes(i),numMatches(k),doExact);
        end
    end
    fprintf(1,'\n');
end
%-------------------------------------------------------------------------------
fileName = fullfile('Data','preComputePositiveMatch.mat');
save(fileName,'probMatches','numOnes','numMatches','whatMax','doExact','numGenes');
fprintf(1,'Saved to %s\n',fileName);

end
