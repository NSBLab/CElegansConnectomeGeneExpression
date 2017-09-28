function [PhiNorm,PhiTrue,PhiRand, rAdj] = RichClubPhiNorm_A(Adj,kmax,numIter,numRepeats,WhatTypeNetwork,whatNullModel, wei_freq)
% RichClubPhiNorm
%
% INPUTS:
% Adjacency matrix, Adj
% Compares to m randomly permuted version of the same matrix
% Does QE iterations, where E is the number of edges
% In the case of weighted networks, THIS WILL NOT preserve strength of each node
% (i.e., weight of its links) but will preserve the degree of each node.
%
% Rewire the network NumLinks*Q times, and repeat this m times.
%
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-16, updated for use in analyzing the mouse connectome.
% Ben Fulcher, 2014-04-28
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs, preliminaries
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(kmax)
    kmax = max(sum(Adj));
    fprintf(1,'Setting maximum k to k_max = %u\n',kmax);
end

if nargin < 6 || isempty(whatNullModel)
    switch WhatTypeNetwork
    case 'bu'
        whatNullModel = 'randmio_und';
    case 'bd'
        whatNullModel = 'randmio_dir';
    case 'wu' 
        whatNullModel = 'Randomise_weights_keeping_topology';
    case 'wus'
        whatNullModel = 'Randomise_weights_keeping_topology';
    otherwise
        error('Appropriate default null model unavailable.');
    end
end

NumNodes = length(Adj);

% ------------------------------------------------------------------------------
% Compute phi for the adjacency matrix
% ------------------------------------------------------------------------------
switch WhatTypeNetwork
case 'wd' % Weighted, directed network:
    RichClubFun = @(x) rich_club_wd(x,kmax);
case 'wu' % Weighted, undirected network:
    RichClubFun = @(x) rich_club_wu(x,kmax);
case 'wus' % Weighted, undirected network:
    RichClubFun = @(x) rich_club_wu_strength(x,kmax);
case 'bu2h' % Weighted, undirected network:
    RichClubFun = @(x) rich_club_bu(x,kmax);
case 'bu' % binary, undirected network:
    RichClubFun = @(x) rich_club_bu(x,kmax);
case 'bd' % binary, directed network:
    % "degree is taken to be the sum of incoming and outgoing connections"
    RichClubFun = @(x) rich_club_bd(x,kmax);
otherwise
    error('Unknown network type ''%s''',WhatTypeNetwork);
end

% ------------------------------------------------------------------------------
% Compute phi for the real network:
% ------------------------------------------------------------------------------
PhiTrue = RichClubFun(Adj);

% ------------------------------------------------------------------------------
% Compute for randomized versions of the network
% ------------------------------------------------------------------------------
PhiRand = zeros(numRepeats,kmax);
fprintf(1,['Computing %u link-permuted matrices, ' ...
                    'using %u iterations for each randomization\n'],numRepeats,numIter);

switch whatNullModel
case 'randmio_und'
    f_rand_null = @randmio_und;
case 'randmio_dir'
    f_rand_null = @randmio_dir;
case 'shuffleWeights' % topology fixed, strength preserved, weight distribution not preserved
    f_rand_null = @f_shuffleWeights;
case 'Randomise_weights_keeping_topology' % topology kept, weights distributed randomly
    f_rand_null = @Randomise_weights_keeping_topology;
case 'null_model_und_sign'
    f_rand_null = @null_model_und_sign;
end

timer = tic;
parfor i = 1:numRepeats
    fprintf(1,'[%u/%u] Rewiring each link %u times...',i,numRepeats,numIter);
    if strcmp(whatNullModel, 'null_model_und_sign')
    [Adj_rand, numRewirings] = f_rand_null(Adj, numIter, wei_freq); % Random graph with preserved in/out degree distribution
    else 
    [Adj_rand, numRewirings] = f_rand_null(Adj, numIter); % Random graph with preserved in/out degree distribution   
    end
    fprintf(1,' %u rewirings performed.\n',numRewirings);
    PhiRand(i,:) = RichClubFun(Adj_rand);

%     if i==1 || mod(i,numRepeats/10)==0
%         fprintf(1,'Approx. %s remaining...\n',BF_thetime(toc(timer)/i*(numRepeats-i)));
%     end
end


% ------------------------------------------------------------------------------
% Calculate normalized phi values for the resulting nomalization, PhiNorm
% ------------------------------------------------------------------------------
% Definition in vandenHeuvel:2011he is the following:
meanNotNaN = @(x)mean(x(~isnan(x) & ~isinf(x)));
PhiNorm = PhiTrue./arrayfun(@(x)meanNotNaN(PhiRand(:,x)),1:kmax);

% ------------------------------------------------------------------------------
% Extra functions:
function [rAdj,numRewirings] = f_shuffleWeights(Adj,numIter);
    % Ben Fulcher, 2014-12-01
    % Shuffles weights, keeping the topology fixed while preserving
    % strenght of each node

    % Get all elements of link data where a connection exists:
%     allActualLinks = Adj(Adj~=0); % Shuffle them 
    [suffled]=random_weight_distribution(Adj);
    
    allActualLinksDataShuffled = suffled;

    % Put them back in the matrix
    rAdj = zeros(size(Adj));
    rAdj(Adj~=0) = allActualLinksDataShuffled;


    % Not relevant to this method:
    numRewirings = 0;
end
% ------------------------------------------------------------------------------

end

