function [PhiNorm,PhiTrue,PhiRand] = RichClubPhiNorm(Adj,kmax,numIter,numRepeats,WhatTypeNetwork,whatNullModel)
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
case 'shuffleWeights'
    f_rand_null = @f_shuffleWeights;
end

timer = tic;
parfor i = 1:numRepeats
    fprintf(1,'[%u/%u] Rewiring each link %u times...',i,numRepeats,numIter);
    [Adj_rand,numRewirings] = f_rand_null(Adj,numIter); % Random graph with preserved in/out degree distribution
    fprintf(1,' %u rewirings performed.\n',numRewirings);
    PhiRand(i,:) = RichClubFun(Adj_rand);

    if i==1 || mod(i,numRepeats/10)==0
        fprintf(1,'Approx. %s remaining...\n',BF_thetime(toc(timer)/i*(numRepeats-i)));
    end
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
    % Shuffles weights, keeping the topology fixed

    % 1. Find
    % Get all elements of link data where a connection exists:
    allActualLinks = Adj(Adj~=0);

    % Shuffle them:
    allActualLinksDataShuffled = allActualLinks(randperm(length(allActualLinks)));

    % Put them back in the matrix
    rAdj = logical(Adj)+0; 
    rAdj(Adj~=0) = allActualLinksDataShuffled;

    % Not relevant to this method:
    numRewirings = 0;
end
% ------------------------------------------------------------------------------

end

% % ------------------------------------------------------------------------------
% % Get normalized phi by comparing to randomized versions
% %    (that preserve the degree distribution)
% % ------------------------------------------------------------------------------
% % Rewire numIter times
% % Do this m times
%
% % Much easier in the format NodeLinks = [Node_i,Node_j,weight];
% % AllLinks = Adj(logical(triu(ones(size(Adj)),+1)));
% % Set below-diagonal terms to zero, assuming a symmetric input matrix
% AllLinks = Adj;
% AllLinks(logical(tril(ones(size(Adj)),-1))) = 0;
% [Thei,Thej] = find(AllLinks>0);
% NumLinks = length(Thei);
% numIter = NumLinks*Q;
%
% NodeLinks = zeros(length(Thei),3);
% for i = 1:NumLinks
%     NodeLinks(i,1) = Thei(i);
%     NodeLinks(i,2) = Thej(i);
%     NodeLinks(i,3) = Adj(Thei(i),Thej(i));
% end
%
% k_orig = (sum(Adj>0)); % Degree of each node
%
% fprintf(1,['Computing %u link-permuted matrices, ' ...
%                     'using %u iterations for each randomization\n'],m,numIter);
% % (as per Julian J. McAuley et al. Appl. Phys. Lett. 91, 084103, 2007)
% % (the swapping method in Milo et al., arXiv:2004)
% PhiRand = zeros(m,kmax);
% for i = 1:m % m random matrices
%     fprintf(1,'\n---------ITERATION%u---------\n',i);
%     NodeLinks_ran = NodeLinks; % start the same, then progressively randomize
%
%     DidSwap = zeros(numIter,1);
%     IterTimer = tic; % time the iterations
%     for j = 1:numIter % numIter link permutations
%
%         if (mod(j,floor(numIter/10))==0)
%             fprintf(1,'Approx %s remaining! We have made %u swaps so far...\n', ...
%                         BF_thetime(toc(IterTimer)/j*(numIter-j)),sum(DidSwap(1:j-1)))
%         end
%
%         % Pick two links at random: r1 and r2
%         r1 = randi(NumLinks);
%         r2 = randi(NumLinks);
%         % % Make sure the second link is not same as the first link:
%         % notr1 = setxor(1:NumLinks,r1);
%         % r2 = notr1(r2);
%
%         % Swapping will either have no effect, or produce a self-link
%         if length(unique([NodeLinks_ran(r1,1:2),NodeLinks_ran(r2,1:2)])) < 4
%             % fprintf(1,'(%u,%u) and (%u,%u) self-link\n', ...
%             %                         NodeLinks_ran(r1,1),NodeLinks_ran(r1,2), ...
%             %                         NodeLinks_ran(r2,1),NodeLinks_ran(r2,2));
%             % fprintf(1,'S');
%             continue
%         end
%
%         % Pick start (or end)-point to swap
%         n1 = randi(2); % swap node (1) or (2) in first link
%         n2 = randi(2); % with node (1) or (2) in the second link
%
%         % Swap node n1 of link r1 with node n2 of link r2
%         old1 = NodeLinks_ran(r1,1:2);
%         old2 = NodeLinks_ran(r2,1:2);
%         new1 = old1; new1(n1) = old2(n2); new1 = sort(new1);
%         new2 = old2; new2(n2) = old1(n1); new2 = sort(new2);
%
%         % Check that these new links don't already exist (else you create a double-link):
%         if any(sum(bsxfun(@minus,NodeLinks_ran(:,1:2),new1)==0,2)==2)
%             % fprintf(1,'Can''t swap (%u,%u) -> (%u,%u) because it already exists...\n', ...
%                                             % old1(1),old1(2),new1(1),new1(2));
%             continue
%         elseif any(sum(bsxfun(@minus,NodeLinks_ran(:,1:2),new2)==0,2)==2)
%             % fprintf(1,'Can''t swap (%u,%u) -> (%u,%u) because it already exists...\n', ...
%                                             % old2(1),old2(2),new2(1),new2(2));
%             % fprintf(1,'E');
%             continue
%         end
%
%         % So we haven't continued yet and can swap!
%         % fprintf(1,'Swapped (%u,%u) -> (%u,%u) and (%u,%u) -> (%u,%u)\n', ...
%         %                         old1(1),old1(2),new1(1),new1(2), ...
%         %                         old2(1),old2(2),new2(1),new2(2));
%         NodeLinks_ran(r1,1:2) = new1;
%         NodeLinks_ran(r2,1:2) = new2;
%         DidSwap(j) = 1;
%
%         % Adj_rand = zeros(NumNodes);
%         % for k = 1:NumLinks
%         %     Adj_rand(NodeLinks_ran(k,1),NodeLinks_ran(k,2)) = NodeLinks_ran(k,3);
%         % end
%         % Adj_rand = Adj_rand + Adj_rand';
%         % k_rand = (sum(Adj_rand>0));
%         % if sum(k_rand-k_orig)~=0
%         %     fprintf(1,'WOWWOWWOW\n');
%         %     keyboard
%         % end
%
%     end
%
%     fprintf(1,'Made a total of %u (/%u) link permutations...\n',sum(DidSwap),numIter);
%
%     % Convert back to an adjacency matrix to feed into the rich club calculation:
%     Adj_rand = zeros(NumNodes);
%     for j = 1:NumLinks
%         Adj_rand(NodeLinks_ran(j,1),NodeLinks_ran(j,2)) = NodeLinks_ran(j,3);
%     end
%     Adj_rand = Adj_rand + Adj_rand';
%     k_rand = (sum(Adj_rand>0));
%     if sum(k_rand-k_orig)~=0
%         fprintf(1,'Sorry to break this to you, but degree distributions actually don''t match... :-/\n');
%         keyboard
%     else
%         fprintf(1,'Degree distributions match!\n');
%     end
%
%     % Visualize the randomized matrix
%     % if i==1
%         % figure('color','w'); box('on');
%         % subplot(2,1,1); pcolor(Adj); shading('flat'); axis square
%         % subplot(2,1,2); pcolor(Adj_rand); shading('flat'); axis square
%     % end
%
%     PhiRand(i,:) = rich_club_wu(Adj_rand,kmax);
% end
