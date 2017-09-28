function [GenesinGO, GOs] = searchGOcategories(C,G, whatToCompare)
%-------------------------------------------------------------------------------
% Aurina Arnatkeviciute 24-04-2017
% This function checks if main GO categories represented in mouse
% (Fulcher,2016) are found in
% 1. worm gene ontology
% 2. among subsample of 948 genes (all)
% 3. among subsample of 577 genes (connected vs unconnected)
% 4. among subsample of 390 genes (RichFeeder vs peripheral)
% also can look for any other type of go categories: 'glutamate',
% 'choline' etc. by tying a string you're interested in in whatTOCompare
% variable.
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
load('GOAnnotation.mat');
load('GOTerms_BP.mat');
% define main GO categories enriched in mouse (Table 1, main text)
switch whatToCompare
    case 'mouse'
        compareGO = [6099;6101;15988;15991; 1902600];
    otherwise
        comp = cell(size(GOTable,1),1);
        for i=1:size(GOTable,1)
            transmiter = sprintf('%s', whatToCompare);
            comp{i} = findstr(transmiter, GOTable.GOName{i});
        end
        overlap = find(~cellfun(@isempty,comp));
        compareGO = GOTable.GOID(overlap);
end

[~, indGO] = intersect(GOTable.GOID, compareGO);
GOnames = GOTable.GOName(indGO);
%-------------------------------------------------------------------------------
%   1. search for those mouse categories in the worm annotation file.
%-------------------------------------------------------------------------------
isinWormGO = intersect(allGOCategories,compareGO);
if isempty(isinWormGO)
    fprintf('NO matching GO categories found\n');
else
%-------------------------------------------------------------------------------
    %2. search for those mouse categories in the subsample of 900 genes
%-------------------------------------------------------------------------------
    %initiate variable
    GOcategories = cell(length(allGOCategories),1); indGOcategories = cell(length(allGOCategories),1);
   % GOs = cell(3,2);
    for k=1:3

        if k==1
            GeneList = G.geneAcronyms.Direct;
        elseif k==2
            [TCON] = GCCbinomial(C,G,'Connected', true, false);
            GeneList = TCON.geneName;
        elseif k==3
            [TRF] = GCCbinomial(C,G,'RichFeeder', true, false);
            GeneList = TRF.geneName;
        end


        for i=1:length(allGOCategories)

            [GOcategories{i}, indGOcategories{i}] = intersect(geneAcronymAnnotations{i}, GeneList);

        end

        inds = find(~cellfun(@isempty,GOcategories));

        GOcats = allGOCategories(inds);
        % get GO names
        [~, GOinds] = intersect(GOTable.GOID, GOcats);
        GOwormID = GOTable.GOID(GOinds);

        % look for intersection

        IsinSubsample = intersect(GOwormID, compareGO);
        [GOs{k,1},indGO]=intersect(allGOCategories,IsinSubsample);
        if isempty(indGO)
            fprintf('NO matching GO categories found\n');
        else
        [~,I]=intersect(GOTable.GOID,GOs{k,1}, 'stable');
        GOs{k,2} = GOTable.GOName(I);
        WhatGenes = geneAcronymAnnotations(indGO);
        end


        % how many of these are in our list?

        for j=1:length(WhatGenes)
            GenesinSubsample{j,k} = intersect(WhatGenes{j},GeneList);
        end
    end
        % put results into the table    
        GenesinGO = table();
        GenesinGO.all948 = GenesinSubsample(:,1);
        GenesinGO.con = GenesinSubsample(:,2);
        GenesinGO.rf = GenesinSubsample(:,3);

end
end
