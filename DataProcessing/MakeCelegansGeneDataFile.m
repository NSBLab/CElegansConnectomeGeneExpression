function G = MakeCelegansGeneDataFile(C, whatData, annotationType)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-10-19
% Creates G structure and saves it to a file according to a dataset and
% annotation type (or a combination of them) chosen
% (WS249, WS254, WS255). Needs C structure and neuron x gene expression
% matrix which is saved in a separate file.
% ------------------------------------------------------------------------------
% This script saves all gene data into one file "CelegansGeneData.mat".
% what data - choose 249 or 254 dataset.
% load connectivity data (will need RegionStructure from this dataset)
% annotationType - insert final annotationType (if a combination is chosen,
% use all of them in the right order without spaces (e.g. CertainEmpty) -
% the way a WSXXX_annotationType.mat file is saved)
%-------------------------------------------------------------------------------

if nargin < 1
    load('CElegansConnectivityData.mat','C');
end
if nargin < 2
    whatData = 256;
end
if nargin < 3
    annotationType = 'CertainEmptyEnrichedPartial';
    fprintf(1,'Using ''CertainEmptyEnrichedPartial'' annotations by default...\n');
end
%-------------------------------------------------------------------------------

% Assign RegionStruct based on connectivity data
G.RegionStruct = C.RegionStruct;

% Load gene data (gene2neuron are assigned using GeneExpressionFinder.m)
fileName = sprintf('WS%d_%s.mat',whatData,annotationType);
load(fileName);
numGenes = size(Direct.expressionMatrix,2);
fprintf(1,'Data from %s loaded (%u genes)...\n',fileName,numGenes);

%[~, iInd]=intersect(Indirect.geneIDsMatrix, Direct.geneIDsMatrix);
% NOTE: number of intersect between indirect and direct is different
% because when indirrect annotations are added, some genes become expressed
% in all neurons and are removed from the final matrix)

% This Direct 0/1 assignment defines if a gene in the Indirect+direct
% matrix was assigned dirrectly or indirrectly.
% we are not using indirectly annotated genes, so only direct data is
% saved.
for l = 1:numGenes
    G.GeneStruct(l).acronym = Direct.geneAcronymMatrix{l};
    G.GeneStruct(l).geneid = Direct.geneIDsMatrix(l);
    G.GeneStruct(l).Direct = true;
end
% assign 1 for directly annotated genes
% for j=1:length(iInd)
% G.GeneStruct(iInd(j),1).Direct=1;
% end

% this number of genes defines directly and indirectly assigned genes
% separately. Direct - only direct. Indirect - (direct+indirect).
G.numGenes.Direct = length(Direct.geneIDsMatrix);
%G.numGenes.Indirect=length(Indirect.geneIDsMatrix);

G.GeneExpData.Direct = Direct.expressionMatrix;
%G.GeneExpData.Indirect=Indirect.expressionMatrix;

G.geneIDs.Direct = Direct.geneIDsMatrix;
%G.geneIDs.Indirect = Indirect.geneIDsMatrix;

G.geneAcronyms.Direct = Direct.geneAcronymMatrix;
%G.geneAcronyms.Indirect = Indirect.geneAcronymMatrix;

end
