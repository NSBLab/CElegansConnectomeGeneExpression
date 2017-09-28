function [expressionMatrix,geneIDsMatrix,geneAcronymMatrix] = GeneExpressionFinder(C,whatData,doDirect,deleteGenes,annotationType1,annotationType2,annotationType3,annotationType4)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-10-19
% this function for each neuron assigns expressed genes.
% ------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% 1.1 Loads connectivity data to structure C
% 1.2 whatData = 249, 254 or 255 (gene data file to use)
% 1.3 doDirect = 1/0; 1 - only direct annotations; 0 - dirrect+undirrect annotations
% 1.4.1. for WS249: annotationType = 'Empty','Certain', 'Life_stage',
% 'Partial', 'Uncertain' 0 - use full dataset
% 1.4.2. for WS254-255: annotationType = 'Empty','Certain', 'Enriched',
% 'Partial', 'Uncertain', 0 - use full dataset
% IMPORTANT:annotationType1-5:
% if annotationType1=0 - all genes; If you want to combine several
% annotation types, write them in any order;

% 2. Loads gene information from a selected file 'anatomy_association.WSXXX.txt'

% 2. Assigns all genes to specific neurons. If a specific neuron is not mentioned under a gene, it goes down hierarchy to find a specific neuron under a parent name.
% 3. Outputs results to matlab file: GeneExpression.mat
%-------------------------------------------------------------------------------

timer = tic;
if nargin < 1
    load('CElegansConnectivityData.mat');
    fprintf(1,'Connectivity loaded\n');
end
if nargin < 2
    whatData = 256;
    fprintf(1,'Using WB256 by default\n');
end
if nargin < 3
    doDirect = true;
    fprintf(1,'Using direct annotations by default\n');
end
if nargin < 4
    deleteGenes = true;
    fprintf(1,'Deleting genes by default\n');
end
if nargin < 5
    annotationType1 = 0;
end
if nargin < 6
    annotationType2 = [];
end
if nargin < 7
    annotationType3 = [];
end
if nargin < 8
    annotationType4 = [];
end

%-------------------------------------------------------------------------------
% Extract neuron ID, count neurons
neuronID = cell2mat({C.RegionStruct.id});
numNeurons = length(neuronID);
fprintf(1,'%u neurons\n',numNeurons);

%-------------------------------------------------------------------------------
%% Read columns of data as strings:
% read in gene file
fileName = sprintf('anatomy_association.WS%d.txt',whatData);

fileID = fopen(fileName,'r'); % read in text file % Formatspec: if an error occurs for a different file, try regenerating the code from the Import Tool.
if fileID==-1
    error('Cannot read in %s (try regenerating the code from the Import Tool)',fileName);
end
fprintf(1,'Loaded anatomy associations from %s\n',fileName);

formatSpec = '%*s%s%s%s%s%s%[^\n\r]';
delimiter = '\t';
% Parse the input text file:
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'ReturnOnError',false);
fclose(fileID);

% Replace IDs with letters with number IDs
f_getIdN = @(x)str2num(x(6:end));
f_getIdG = @(x)str2num(x(7:end));

fprintf('Changing anatomy term and gene IDs to numerical\n')
dataArrayMod{1} = cellfun(f_getIdG,dataArray{1});
geneIDs = dataArrayMod{1};
dataArrayMod{3} = cellfun(f_getIdN,dataArray{4});
anatomyTerm = dataArrayMod{3};
dataArrayMod{2} = dataArray{2};
geneAcronym = dataArrayMod{2};
%geneNames = dataArrayMod{2};
expressionStatus = dataArray{3};
% replace all blank cells with a string "Empty" - makes further filtering
% easier.
emptyStatus = find(cellfun(@isempty,expressionStatus));
expressionStatus(emptyStatus) = {'Empty'};

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% get indexes for a certain type of status "certain","uncertain",'enriched'
% (or combination of them)
fprintf('Filtering genes according to expression status annotation\n' )

if sum(annotationType1~=0)~=0
    compareAnnotations(:,1) = strcmp(annotationType1,expressionStatus);
    compareAnnotations(:,2) = strcmp(annotationType2,expressionStatus);
    compareAnnotations(:,3) = strcmp(annotationType3,expressionStatus);
    compareAnnotations(:,4) = strcmp(annotationType4,expressionStatus);

    sumAnnot = sum(compareAnnotations,2);
    INDstatus = find(sumAnnot==1);
    geneIDs = geneIDs(INDstatus);
    anatomyTerm = anatomyTerm(INDstatus);
    geneAcronym = geneAcronym(INDstatus);
end

uniqueGenesID = unique(geneIDs, 'stable');
uniqueGenesAcronym = unique(geneAcronym, 'stable');
numGenes = length(uniqueGenesID);

fprintf(1,'%u genes annotated \n',numGenes);
%-------------------------------------------------------------------------------
%% Create lists for hierarchy - parent->child
%-------------------------------------------------------------------------------
% Read in and parse data from hierarchy.csv
fprintf('Reading hierarchy file\n');
[parentName,childName,parentID,childID] = ReadInHierarchy();

%-------------------------------------------------------------------------------
%% For each gene find a neuron it was DIRECTLY (or DIRECTLY+INDIRECTLY) annotated to
%% Remove genes that were not annotated or were annotated to all neurons
%-------------------------------------------------------------------------------
expressionMatrix = fillMatrix(); %function below
fprintf('Expression matrix created\n')
if isempty(expressionMatrix)
    fprintf('No genes with this annotation were assigned\nc')
else
    % remove genes that are not expressed in any neurons or are expressed in
    % all neurons if several annotations are chosen. If one annotation is
    % chosen, keep all entries in case you'll need to combine data afterwards,
    fprintf('Filterring out genes that expressed in all or none neurons\n')
    geneIDsMatrix = uniqueGenesID;
    geneAcronymMatrix = uniqueGenesAcronym;
    if deleteGenes == 1
        sumMatrix = sum(expressionMatrix,1);
        zeroInd = find(sumMatrix == 0);
        allInd = find(sumMatrix == 279);
        indExclude = horzcat(zeroInd,allInd);

        expressionMatrix(:, indExclude) = [];
        geneIDsMatrix(indExclude) = [];
        geneAcronymMatrix(indExclude) = [];
    end
end

fprintf(1,'neurons assigned to genes in %s\n',BF_thetime(toc(timer)));

%-------------------------------------------------------------------------------
    function expressionMatrix = fillMatrix()
        expressionMatrix = zeros(numNeurons,numGenes);

        for m = 1:numGenes
            INDAnatomyGene = geneIDs==uniqueGenesID(m);
            subsetForGene = anatomyTerm(INDAnatomyGene);

            neuronsForGene = zeros(numNeurons,1);
            for l=1:length(subsetForGene)
                neuronsGeneL = PropagateDown(subsetForGene(l));
                neuronsForGene(neuronsGeneL==1) = 1;
            end
            expressionMatrix(:,m) = neuronsForGene;
        end
    end

    function neurons = PropagateDown(idStart)
        % start at an ID and propagate down the hierarchy, annotating matches on the way
        idHere = idStart;
        neurons = zeros(numNeurons,1);

        % find matches to current ID
        if doDirect==1
            if any(idHere==neuronID)
                neurons(idHere==neuronID) = 1;
            %elseif any(idHere==Parents)
                %neurons(idHere==Parents) = 1;
            end
        else
            if any(idHere==neuronID)
                neurons(idHere==neuronID) = 1;
                %fprintf(1,'Annotated for %u\n',idStart);
            else
                % Go down hierarchy to children
                idDown = childID(parentID==idHere);
                numChildren = length(idDown);
                if numChildren == 0
                    %fprintf(1,'Got to bottom for %u\n',idStart);
                else
                    %fprintf(1,'Going down for %u children\n',numChildren);
                    for j = 1:numChildren
                        % go down hierarchy for each parent:
                    neuron_j = PropagateDown(idDown(j));
                    end
                    neurons(neuron_j==1) = 1;
                end
            end
        end
    end

end
