% get expression data from O.Hobert's paper

clear all; close all; 
cd('/Users/Aurina/GoogleDrive/Genetics_connectome/CElegans/CElegansCode/Data');
load('CElegansConnectivityData.mat'); 
DATA = textread('mmc2mod.txt', '%s','delimiter', '\n');
for ii=1:length(DATA)
    N = DATA{ii,1};
    [neuron, genes] = strtok(N); 
    Neurons{ii} = neuron; 
    Genes{ii,:} = genes; 
end
% order of neurons is the same, so just use neuron names from C. 
 Neurons = C.RegionAcronyms; 
 
 for jj=1:279
     genelist = Genes{jj,:}; 
     kk=1; 
     while ~isempty(genelist)
     
         [selectedGene, genes] = strtok(genelist); 
         O{jj,kk} = selectedGene; 
         genelist = genes;
         kk=kk+1; 
     end
 end
 
 On = O(:,1:141);% delete all empty columns 
 R = On(:); % arrange all into one column
 Ugenes = uniquecell(R); % get unique genes; 
 Ugenes = Ugenes(~cellfun('isempty',Ugenes)); 
 
 % arragne data into a matrix
 geneExpression = zeros(279,931); 
 for i = 1:length(Ugenes)
    for j = 1:279
        if any(strcmp(On(j,:),Ugenes{i}))
            geneExpression(j,i) = 1;
        end
    end
 end

% aragne data into G structure; 
G.RegionStruct = C.RegionStruct; 
G.numGenes.Direct = length(Ugenes);
G.GeneExpData.Direct = geneExpression; 
G.geneAcronyms.Direct = Ugenes;

save('HobertGeneData.mat','G'); 
    