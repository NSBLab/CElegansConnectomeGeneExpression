% check if there are genes in worm from the categories used in mouse
% get a list of go categories in mouse and compare to the list in worm (for
% all posible genes)
% ------------------------------------------------------------------------------
% read in all txt files from mouse and combine unique GO categories in a
% list
function GOtable = compareGOannotations(C,G)
%-------------------------------------------------------------------------------
% Aurina Arnatkeviciute 13-04-2017
%
% This function compares GO annotations found in mouse by Fulcher PNAS 2016 to GO annotations available in C. elegans"
% 1. to all possible GO annotations in C elegans (from GO annotation file for the worm - all possible annotations generated using
% 2. to GO annotations for genes in C.elegans dataset (948 genes)
% 3. to GO annotations for genes in C.elegans dataset that were used for connected vs unconnected analysis (577 genes)
% 4. to GO annotations for genes in C.elegans dataset that were used for rich/feeder vs peripheral analysis (390 genes)
% INPUT: C and G structures.
% OUTPUT: a list of genes for each GO category that was matched geneListinSample (948, 577,390 gene samples used)
%-------------------------------------------------------------------------------

load('GOAnnotation.mat');
types = {'conBP', 'conBPCC', 'rfBP', 'rfBPCC'};
GOmouse = cell(length(types),1);
for t=1:length(types)
    filePath = sprintf('Data/ermineJdata/mouseGO%s.txt', types{t});
    fid = fopen(filePath,'r');
    % headings = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s',1,'Delimiter','\t','CollectOutput',1);
    GOmouse{t} = textscan(fid,'%s%s%*[^\n\r]','Delimiter','\t'); %,'EndOfLine','\n');
    fclose(fid);
    %textscan(fid,'%s%s%s%*[^\n\r]','Delimiter','  & ','EndOfLine','\n');
end

GOids = [GOmouse{1}{1,1};GOmouse{2}{1,1};GOmouse{3}{1,1};GOmouse{4}{1,1}];
GOnames = [GOmouse{1}{1,2};GOmouse{2}{1,2};GOmouse{3}{1,2};GOmouse{4}{1,2}];

f_ind = @(x)str2num(x(5:end));
GOlistmouse = cellfun(f_ind,GOids);

[~, indGOmouseANDworm] = intersect(allGOCategories, GOlistmouse);
[~, indS] = intersect(GOlistmouse, allGOCategories);
GOmouseANDwormNames = GOnames(indS);

% ------------------------------------------------------------------------------
% get genes with matching categories in worm
genesWORM = geneAcronymAnnotations(indGOmouseANDworm);
% for each GO category, check how many genes we have in our list of 948
% genes.
% all genes
geneListALL = G.geneAcronyms.Direct;
% genes for connected vs unconnected analysis
[~,geneScoresCON,~] = GCCbinomial(C,G,'Connected',true,false);
geneListCON= geneScoresCON.geneName;
%genes for Rich feeder analysis
[~,geneScoresRF,~] = GCCbinomial(C,G,'RichFeeder',true,false);
geneListRF= geneScoresRF.geneName;

geneListinSample = cell(length(genesWORM),4);
analysis = {'all', 'connected', 'richfeeder'};
for anal=analysis
    for j=1:length(genesWORM)
        GOcategory = genesWORM{j};
        if strcmp(anal,'all')
            geneListinSample{j,1} = intersect(GOcategory, geneListALL);
        elseif strcmp(anal,'connected')
            geneListinSample{j,2} = intersect(GOcategory, geneListCON);
        elseif strcmp(anal,'richfeeder')
            geneListinSample{j,3} = intersect(GOcategory, geneListRF);
        end
    end
end
geneListinSample(:,4) = GOmouseANDwormNames;

GOtable = table();
GOtable.list948genesALL = geneListinSample(:,1);
GOtable.list577genesCON = geneListinSample(:,2);
GOtable.list390genesRF = geneListinSample(:,3);
GOtable.GOcategories = geneListinSample(:,4);


filePath = fullfile('Data','ermineJdata','GOannotationsCompared.mat');
save(filePath,'GOtable');
fprintf(1,'Saved to %s\n',filePath);

 %-------------------------------------------------------------------------------
 %-------------------------------------------------------------------------------


end
