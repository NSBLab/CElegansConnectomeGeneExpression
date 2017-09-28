function C = Print_neuronconnect()
%% Aurina 22/09/2016
% In this script several scripts are brought together in order to gather information about Celegans connectivity data from the original .csv file to matlab format.
% Some corrections for easier file manipulation are added to original scripts.

%% Notes
% In a case of chemical synapses, only R+Rp connections are considered (S+Sp = R+Rp): N1 receives synapses from N2, it means that N2 is the source and N1 is the target.
% In order to avoid another step of reordering columns, we write them into a txt. file in this order: N2 (source), N1 (target), Nb(weight) for chemical synapses.

%-------------------------------------------------------------------------------
%% Step 1: Create separate lists of connections for chemical and electrical synapses and save them in separate .txt files
%-------------------------------------------------------------------------------
% go to where the original .csv file is located

C = struct();
inputFileName = 'NeuronConnect(corrected).csv';
fprintf(1,'Reading in connectivity data from %s\n',which(inputFileName));
fid = fopen(inputFileName);
D = textscan(fid, '%s%s%s%s','Delimiter',',');
fclose(fid);

N1 = D{1}; % neuron (1) (pre-synaptic)
N2 = D{2}; % neuron (2) (post-synaptic)
Tp = D{3};
Nb = D{4};

%===============================================================================
%% Create separate lists of connections for chemical and electrical synapses
%===============================================================================

chemicalData = fullfile('Data','Chemical_Celegans279.txt');
electricalData = fullfile('Data','EJ_Celegans279.txt');

% Prints only the R+Rp (=S+Sp) connections: N1 receives synapses from N2, therefore source is N2 and target is N1
fid = fopen(chemicalData,'w');
for i = 1:length(N1)
    if strcmp(Tp(i),'R') || strcmp(Tp(i),'Rp')
        fprintf(fid,'%s\t%s\t%s\t\n',N2{i},N1{i},Nb{i});
    end
end
fclose(fid);

% Prints only the Ej connections (symetrical - do difference which one is source and target)
fid = fopen(electricalData,'w');
for i = 1:length(N1)
    if strcmp(Tp(i),'EJ')
        fprintf(fid,'%s\t%s\t%s\t\n',N1{i},N2{i},Nb{i});
    end
end
fclose(fid);

%-------------------------------------------------------------------------------
%% Arrange connectivity data to matrices for chemical, electrical and chemical+electrical connections.
%-------------------------------------------------------------------------------

% Neuron 1 - source, Neuron 2 - traget.
[Neuron1, Neuron2, weight] = textread(chemicalData,'%s %s %s');

% list of region names:
load('Celegans_positions.mat');

unique_regions_t = Worm279_labels;
unique_regions_s = Worm279_labels;

numNeurons = length(unique_regions_s);

fprintf(1,'We have %u unique neurons!\n',numNeurons)

% Prepare matrices for all connectomes.
Adj_wei_Ch = zeros(length(unique_regions_s),length(unique_regions_t));
Adj_bin_Ch = zeros(length(unique_regions_s),length(unique_regions_t));
Adj_wei_EJ = zeros(length(unique_regions_s),length(unique_regions_t));
Adj_bin_EJ = zeros(length(unique_regions_s),length(unique_regions_t));
Adj_wei_all = zeros(length(unique_regions_s),length(unique_regions_t));
Adj_bin_all = zeros(length(unique_regions_s),length(unique_regions_t));

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%% CHEMICAL:
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% weighted (R+Rp)
fprintf('Arranging chemical weighted matrix\n')
k=1;
for i = 1:length(unique_regions_s)
    for j = 1:length(unique_regions_t)
        if (strcmp(Neuron1{k},unique_regions_s{i}) && strcmp(Neuron2{k},unique_regions_t{j}))
            if k~=length(Neuron1)
                if (strcmp(Neuron1{k}, Neuron1{k+1}) && strcmp(Neuron2{k}, Neuron2{k+1}))
                    Adj_wei_Ch(i,j) = str2double(weight{k})+str2double(weight{k+1}); %R + Rp
                    if k < length(Neuron1)-1
                        k = k+2;
                    end
                else
                    Adj_wei_Ch(i,j) = str2double(weight{k});
                    if k < length(Neuron1)
                        k = k+1;
                    end
                end
            else
                Adj_wei_Ch(i,j) = str2double(weight{k});
                if k < length(Neuron1)
                    k = k+1;
                end
            end
        end
    end
end
Adj_wei_Ch(eye(size(Adj_wei_Ch))~=0)=0;
% unweighted
fprintf('Arranging chemical unweighted matrix\n')
for i = 1:length(unique_regions_s)
    for j = 1:length(unique_regions_t)
        if any(strcmp(Neuron1,unique_regions_s{i}) & strcmp(Neuron2,unique_regions_t{j}))
            Adj_bin_Ch(i,j) = 1;
        end
    end
end
Adj_bin_Ch(eye(size(Adj_bin_Ch))~=0)=0;

%-------------------------------------------------------------------------------
%% ELECTRICAL
%-------------------------------------------------------------------------------
% Neuron 1 - source, Neuron 2 - target.
[Neuron1, Neuron2, weight] = textread(electricalData,'%s %s %s');

fprintf('Arranging electrical weighted and unweighted matrices\n')
for i = 1:length(unique_regions_s)
    for j = 1:length(unique_regions_t)
        for k = 1:length(weight)
            if (strcmp(Neuron1{k},unique_regions_s{i}) && strcmp(Neuron2{k},unique_regions_t{j}))
                % weighted
                Adj_wei_EJ(i,j) = str2double(weight{k});
                % unweighted
                Adj_bin_EJ(i,j) = 1;
            end
        end
    end
end
Adj_bin_EJ(eye(size(Adj_bin_EJ))~=0)=0;
Adj_wei_EJ(eye(size(Adj_wei_EJ))~=0)=0;

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%% CHEMICAL + ELECTRICAL
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% unWeighted
fprintf('Arranging chemical+electrical weighted and unweighted matrices\n')
for i = 1:length(unique_regions_s)
    for j = 1:length(unique_regions_s)
        if (Adj_bin_EJ(i,j)==1) || (Adj_bin_Ch(i,j)==1)
            Adj_bin_all(i,j) = 1;
        end
    end
end

% weighted
for i = 1:length(unique_regions_s)
    for j = 1:length(unique_regions_s)
        Adj_wei_all(i,j) = Adj_wei_EJ(i,j)+Adj_wei_Ch(i,j);
    end
end



%% Step 3: save all matrices in separate .mat files
% unweighted
% save(fullfile('Data','Adj_bin_all.mat'), 'Adj_bin_all');
% save(fullfile('Data','Adj_bin_Ch.mat'), 'Adj_bin_Ch');
% save(fullfile('Data','Adj_bin_EJ.mat'), 'Adj_bin_EJ');
% % weighted
% save(fullfile('Data','Adj_wei_all.mat'), 'Adj_wei_all');
% save(fullfile('Data','Adj_wei_Ch.mat'), 'Adj_wei_Ch');
% save(fullfile('Data','Adj_wei_EJ.mat'), 'Adj_wei_EJ');

%% Step 3: save relevant data into C structure

% matrices
C.Adj_B{1,1} = Adj_bin_Ch;
C.Adj_B_Names{1,1} = 'Adj_bin_Ch';
C.Adj_B{1,2} = Adj_bin_EJ;
C.Adj_B_Names{1,2} = 'Adj_bin_EJ';
C.Adj_B{1,3} = Adj_bin_all;
C.Adj_B_Names{1,3} = 'Adj_bin_all';

C.Adj_W{1,1} = Adj_wei_Ch;
C.Adj_W_Names{1,1} = 'Adj_wei_Ch';
C.Adj_W{1,2} = Adj_wei_EJ;
C.Adj_W_Names{1,2} = 'Adj_wei_EJ';
C.Adj_W{1,3} = Adj_wei_all;
C.Adj_W_Names{1,3} = 'Adj_wei_all';

% Two-dimensional coordinates
C.RegionAcronyms = unique_regions_s;
C.numNeurons = numNeurons;
C.Pos = Worm279_positions;

% label left-side neurons
isLeft = false(279,1);
regionIndex = zeros(279,1);
for n=1:279
    name = C.RegionAcronyms{n};
    regionIndex(n) = n;
    if strcmp(name(end),'L')
        isLeft(n) = true;
    end
end

NeuronNames = unique(Worm279_labels);

% reading hierarchy file
[parentName,childName,parentID,childID] = ReadInHierarchy();

% find IDs for all neurons in neuron list
IDs = vertcat(parentID,childID);
Names = vertcat(parentName,childName);

[Neurons,NeuronIndex] = intersect(Names(:),NeuronNames);
C.NeuronNames = unique(Worm279_labels);
Neuronlist279IDs = IDs(NeuronIndex,:);

% put information to structure
for l = 1:length(unique_regions_s)
    C.RegionStruct(l).acronym = Neurons{l};
    C.RegionStruct(l).id = Neuronlist279IDs(l);
    C.RegionStruct(l).isLeft = isLeft(l);
    C.RegionStruct(l).regionIndex = regionIndex(l);
end

C.Neuronlist279 = [Neurons, num2cell(Neuronlist279IDs)];

%===============================================================================
%% Add lineage data
%===============================================================================
fprintf('Arranging lineage data\n')

fileName1 = 'NeuronLineage_Part1.txt';
fileName2 = 'NeuronLineage_Part2.txt';

D1 = tdfread(fileName1,'\t');
D2 = tdfread(fileName2,'\t');

% combine data from two files
N1 = [cellstr(D1.Neuron_1); cellstr(D2.Neuron_1)]; % neuron 1
N2 = [cellstr(D1.Neuron_2); cellstr(D2.Neuron_2)]; % neuron 2
R = [D1.Relatedness; D2.Relatedness]; % relatedness of neurons 1 and 2

Neurons = C.RegionAcronyms;
LineageDistance = zeros(C.numNeurons, C.numNeurons);
% Assign relatedness data for each pair of neurons.
numAnnotations = length(N1);
for i = 1:numAnnotations
    ix_i = strcmp(Neurons,N1{i});
    ix_j = strcmp(Neurons,N2{i});
    LineageDistance(ix_i,ix_j) = R(i);
end
% LineageDistance = LineageDistance + LineageDistance'; % add another half
% LineageDistance(logical(eye(size(LineageDistance))))=0;
C.LineageDistance = LineageDistance;  % add to structure
% add neuron birth time information

fprintf('Arranging neuron birth time data\n')
load('celegans279_BT.mat');
Names = C.RegionAcronyms;
[~, ind] = intersect(celegans279labels, Names);
% reorder neurons according to standard
birthTimes = celegans279_birth_time(ind);

birthTimeDiff = zeros(279,279);
for i=1:length(unique_regions_s)
    for j=1:length(unique_regions_s)
        birthTimeDiff(i,j) = abs(birthTimes(j)-birthTimes(i));
    end
end

C.BirthTimeDiff = birthTimeDiff;
C.BirthTime = birthTimes;

% make mask for removing distance, lineage and birth time difference information for L/R pairs of neurons
LRmask = GiveMeLRMask(C);
%-------------------------------------------------------------------------------
% Add distance matrix for 2D and 3D coordinates
%-------------------------------------------------------------------------------
%% calculate euclidean distance between neurons
C.Eucl_Dist_2D = squareform(pdist(Worm279_positions,'Euclidean'));
C.Pos3D = CoordinateTo3d(Worm279_positions,isLeft);
C.Eucl_Dist_3D = squareform(pdist(C.Pos3D,'Euclidean'));

C.LineageDistance_noLR = LineageDistance;
C.LineageDistance_noLR(LRmask) = NaN;
C.Eucl_Dist_2D_noLR = C.Eucl_Dist_2D;
C.Eucl_Dist_2D_noLR(LRmask) = NaN;
C.Eucl_Dist_3D_noLR = C.Eucl_Dist_3D;
C.Eucl_Dist_3D_noLR(LRmask) = NaN;
C.BirthTimeDiff_noLR = birthTimeDiff;
C.BirthTimeDiff_noLR(LRmask) = NaN;

%-------------------------------------------------------------------------------
% Add hierarchical information
%-------------------------------------------------------------------------------
fprintf(1,'Adding hierarchical information...');
C = NeuronRegionFinder(C,false);
fprintf(1,' Done.\n');

%-------------------------------------------------------------------------------
% SAVE file
%-------------------------------------------------------------------------------
connDataFile = fullfile('Data','CElegansConnectivityData.mat');
save(connDataFile,'C');
fprintf(1,'Saved connectivity data and annotations to %s\n',connDataFile);
