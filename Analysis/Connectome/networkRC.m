function networkRC(C,synapseType,networkType,nullModel,numIter,numNullNetworks)
% ------------------------------------------------------------------------------
% function calculates and plots normalised rich club coefficient as a
% function of degree.
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C - connectivity structure
% synapse type - chemical 'ch', 'el', 'all';
% networkType - 'bd' - binary directed, 'wd' - weighted directed, 'wu' - weighted undirected, 'bu' - binary undifected.
% nullModel - 'randmio_dir' (chuffle topology), 'shuffleWeights' - shuffle weights, keeping topology.
%   randmio_dir on weighted matrix will provide mixed RC, randmio_dir on
%   binary - will provide topological RC.
% numIter - how many times a link is rewired.                   default - 50
% numNullNetworks - how many null networks is required          default - 1000
% ------------------------------------------------------------------------------
% Outputs
%-------------------------------------------------------------------------------
% <<SAVES TO .mat file in the Data directory>>
% ------------------------------------------------------------------------------
% Normalised RC coefiecient plot
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2017-03-20
% ------------------------------------------------------------------------------

%===============================================================================
% Set defaults
%===============================================================================

if nargin < 2
    synapseType = 'all'; % 'all' - chemical + electrical, 'el' - electrical
end

if nargin < 3
    networkType = 'bd'; % wd - weighted directed, bu -
    fprintf('Binary directed network by DEFAULT\n')
end

if nargin < 4
    nullModel = 'randmio_dir';
    fprintf('Topological RC by DEFAULT\n')
    % whatNullModel = 'shuffleWeights'; - Shuffles weights, keeping the topology fixed
    % whatNullModel = 'randmio_dir'; preserves degree and weight distribution; changes topology. (to be used on binary networks)
end

if nargin < 5
    numIter = 50;
    fprintf(sprintf('%d iterations by DEFAULT\n',numIter))
end

if nargin < 6
    numNullNetworks = 100;
    fprintf(sprintf('%d null networks by DEFAULT\n',numNullNetworks))
end
%===============================================================================
%===============================================================================


Adj = GiveMeAdj(C,'zeroWeighted',synapseType);
[~,~, deg] = degrees_dir(Adj);
kmax = max(deg);

if strcmp(networkType, 'bu') || strcmp(networkType, 'bd' )
    Adj = logical(Adj);
end

% selecting the right null model
if strcmp(synapseType, 'ch')  && strcmp(nullModel, 'randmio_und')
     nullModel = 'randmio_dir';
     fprintf(2,'WRONG null model chosen, changing to directed\n');
elseif strcmp(synapseType, 'all')  && strcmp(nullModel, 'randmio_und')
     nullModel = 'randmio_dir';
     fprintf(2,'WRONG null model chosen, changing to directed\n');
end


% selecting the right network type
if strcmp(synapseType, 'ch')  && strcmp(networkType, 'wu')
    networkType = 'wd';
    fprintf(2,'WRONG network type chosen, changing to directed\n');
elseif strcmp(synapseType, 'all')  && strcmp(networkType, 'wu')
    networkType = 'wd';
    fprintf(2,'WRONG network type chosen, changing to directed\n');
elseif strcmp(synapseType, 'ch')  && strcmp(networkType, 'bu')
    networkType = 'bd';
    fprintf(2,'WRONG network type chosen, changing to directed\n');
elseif strcmp(synapseType, 'all')  && strcmp(networkType, 'bu')
    networkType = 'bd';
    fprintf(2,'WRONG network type chosen, changing to directed\n');
end

% selecting the right null model according to network type
if strcmp(networkType, 'bu') && strcmp(nullModel, 'randmio_dir')
    nullModel = 'randmio_und';
    fprintf(2,'WRONG null model chosen, changing to undirected\n');
elseif strcmp(networkType, 'wu') && strcmp(nullModel, 'randmio_dir')
    nullModel = 'randmio_und';
    fprintf(2,'WRONG null model chosen, changing to undirected\n');
end


if strcmp(networkType, 'bu') || strcmp(networkType, 'bd' )
    Adj = logical(Adj);
end

%-------------------------------------------------------------------------------
% Perform RC coefficient calculation
%-------------------------------------------------------------------------------

[~,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax,numIter,numNullNetworks,networkType,nullModel);

%-------------------------------------------------------------------------------
% Compute p values
%-------------------------------------------------------------------------------
pValues = zeros(kmax,1);
for i = 1:kmax
    pValues(i) = mean(PhiTrue(i) <= PhiRand(:,i));
end

isSig = (pValues <= 0.05);
PhiNormMean = zeros(size(PhiTrue));
for i = 1:size(PhiTrue,2)
    PhiNormMean(i) = PhiTrue(i)/nanmean(PhiRand(:,i));
end

%-------------------------------------------------------------------------------
% Save to file
%-------------------------------------------------------------------------------
fileNameOut = fullfile('Data',sprintf('richClubResults_%s-%s-%s_%u-%uiter.mat',synapseType,...
                                            networkType,nullModel, numNullNetworks,numIter));
save(fileNameOut);
fprintf(1,'Saved results to %s\n',fileNameOut);

end
