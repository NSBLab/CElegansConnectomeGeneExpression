%-------------------------------------------------------------------------------
%% get relatedness matrix (lineage distance) and save to C structure in file
%-------------------------------------------------------------------------------
% load lineage data from files (files saves as txt files instead of the
% original version of xls files)

fileName1 = 'NeuronLineage_Part1.txt';
fileName2 = 'NeuronLineage_Part2.txt';
D1 = tdfread(fileName1,'\t');
D2 = tdfread(fileName2,'\t');

% combine data from two files
N1 = [cellstr(D1.Neuron_1); cellstr(D2.Neuron_1)] ;
N2 = [cellstr(D1.Neuron_2); cellstr(D2.Neuron_2)] ;
R = [D1.Relatedness; D2.Relatedness];

% load connectivity data
fileName = 'CElegansConnectivityData.mat';
load(fileName,'C');

Neurons = C.RegionAcronyms;
LineageDistance = zeros(C.numNeurons, C.numNeurons);

% assign relatedness data for each pair of neurons.
for i=1:C.numNeurons
    for j=1:C.numNeurons
        for k=1:length(N1)
            if (strcmp(N1{k},Neurons{i}) && strcmp(N2{k},Neurons{j}))
                LineageDistance(i,j) = R(k);
            end
        end
    end
end
C.relatedness = LineageDistance;
save(fileName,'C','-append');
