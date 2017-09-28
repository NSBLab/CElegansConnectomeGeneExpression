function runCorrelation(ID)

rng(ID);


minProp = 1; % min number of ones in a vector
maxProp = 150; % max number of ones in a vector - proportion for our data is 0.144
numGenes = 948; % length(G.GeneStruct);
numValues = 150; % number of vectors to calculate coexpression on.

numCases = round(linspace(minProp,maxProp,numValues));

M1 = zeros(length(numCases),numGenes);
for j=1:length(numCases)
    m = numCases(j);
    randVect = zeros(numGenes,1);
    randVect(randperm(numGenes,m)) = 1;
    M1(j,:) = randVect;
end

% calculate correlations between those vectors
numVectors = length(numCases);
BinCorrelation = zeros(numVectors,numVectors);

M2 = squeeze(M1(:,:));
BinCorrelation(:,:) = CoExpressionPositiveMatch(M2,true);
filename = sprintf('%d.mat', ID);
cd ('/projects/kg98/aurina/correlation/');
save(filename, 'BinCorrelation');

end
