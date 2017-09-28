function G = SaveExpressionData(C,whatData)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-10-19
% assigns genes to neurons and saves data for direct
% annotations in a file.
% ------------------------------------------------------------------------------
% choose what gene dataset to use
% whatData = 256, 254, 255; % choose option for the wormbase version.
% choose annotation type (Certain, Partial, Uncertain, Empty, 0 - all
% genes)
%-------------------------------------------------------------------------------

if nargin < 1
    load('CElegansConnectivityData.mat','C');
    fprintf(1,'Loaded connectivity info from ''CElegansConnectivityData.mat''\n');
end
if nargin < 2
    whatData = 256;
    fprintf(1,'Using WB256 BY DEFAULT\n');
end
%-------------------------------------------------------------------------------

annotationType1 = 'Certain';
annotationType2 = 'Empty';
annotationType3 = 'Enriched';
annotationType4 = 'Partial';

% run GeneExpressionFinder function for direct annotations
Direct = struct();
[Direct.expressionMatrix,Direct.geneIDsMatrix,Direct.geneAcronymMatrix] = ...
                GeneExpressionFinder(C,whatData,1,1,annotationType1,annotationType2,annotationType3,annotationType4);

%-------------------------------------------------------------------------------
% Save in data directory
%-------------------------------------------------------------------------------
fileName = sprintf('WS%d_%s%s%s%s.mat',whatData,annotationType1,annotationType2,annotationType3,annotationType4);
fileName = fullfile('Data',fileName);
save(fileName,'Direct');
fprintf(1,'Saved direct annotations to %s\n',fileName);

% uses MakeCelegansGeneDataFile function and creates CelegansGeneData...mat files for each dataset and annotation
% combination
annotationType = strcat(annotationType1,annotationType2,annotationType3,annotationType4);
G = MakeCelegansGeneDataFile(C, whatData, annotationType);
G = Gene_Correlation(G, whatData, annotationType);
G = GeneName2EntrezID(G, whatData, annotationType);

end
