% ------------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Add paths
%-------------------------------------------------------------------------------
fprintf(1,'Adding all subdirectories to the Matlab path...');
% Add paths required for the project (ignoring hidden, including version control)
files = dir;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
if isfield(directories,'folder')
    paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);
else
    paths = arrayfun(@(x)fullfile(pwd,directories(x).name),1:length(directories),'UniformOutput',false);
end
for j = 1:length(paths)
    addpath(genpath(paths{j}))
end
fprintf(1,' Added.\n');

%-------------------------------------------------------------------------------
%% Load data: connectivity matrix and gene expression data
%-------------------------------------------------------------------------------
% cd('/Users/Aurina/GoogleDrive/Genetics_connectome/CElegans/CElegansCode/Data');
fprintf(1,'Loading structural connectivity and gene expression data...');
load('CElegansConnectivityData.mat','C'); % C stands for connectome data
load('CelegansGeneDataWS256CertainEmptyEnrichedPartial.mat','G'); % G stands for gene expression data
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% RegionStruct should match between the G and C datasets
% ------------------------------------------------------------------------------
RegionStruct = G.RegionStruct;

if ~all(arrayfun(@(x)strcmp(C.RegionAcronyms(x),RegionStruct(x).acronym),1:length(C.RegionAcronyms)))
    % Some mismatches between regions in C and G datasets
    error('Regions do not match... :(');
end
