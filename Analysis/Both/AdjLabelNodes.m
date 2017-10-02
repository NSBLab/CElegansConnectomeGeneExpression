% Ben Fulcher, 2014-10-01
% Labels nodes according to their links to/from the rich club
% Modified by Taqi Ali 17-7-15
% ------------------------------------------------------------------------------

function [nodeLabels,labelNames,nodeData] = AdjLabelNodes(labelWhat,Adj,extraParam,networkType)

if nargin < 4
    networkType = 'bd';
end

numNodes = length(Adj);

% ------------------------------------------------------------------------------
% Basic functions:
% ------------------------------------------------------------------------------
% Total degree:

switch networkType
    case {'bu','wu'}
        totDegree = @(x) sum(x,2); % just choose one sided connections
    case {'bd','wd'}
        totDegree = @(x) (sum(x,1)' + sum(x,2)); % outputs a column vector
end


% ------------------------------------------------------------------------------
% Do the labeling:
% ------------------------------------------------------------------------------
switch labelWhat
case {'hub-kth','hub-topN'}
    if strcmp(labelWhat,'hub-topN')
        warning('hub-topN doesn''t work the best yet');
    end
    if nargin < 3
        extraParam = {'degree',50};
    end
    
    whatDegree = extraParam{1};
    switch whatDegree
    case 'degree'
        nodeData = totDegree(Adj);
    end
    
    kth = extraParam{2};
    
    isHub = (nodeData > kth);
    
    
    nodeLabels = isHub + 1;
    
    labelNames = {'non-hub',sprintf('hub (k > %u)',kth)};
    
case 'RichFeed'
    % Ben Fulcher, 2014-10-01
    kth = extraParam{2}; % threshold on degree to count as rich

    % 1. Get total degree of each node:
    nodeData = totDegree(Adj);
    isHub = (nodeData >= kth);

    nodeLabels = zeros(numNodes,1);

    linksToHub = (Adj*isHub > 0);
    linksFromHub = (isHub'*Adj > 0)';


    nodeLabels(isHub) = 1; % Hubs
    nodeLabels(~isHub & linksToHub & ~linksFromHub) = 2; % Feed-in
    nodeLabels(~isHub & ~linksToHub & linksFromHub) = 3; % Feed-out
    nodeLabels(~isHub & linksToHub & linksFromHub) = 4; % Feed-in-out
    nodeLabels(~isHub & ~linksToHub & ~linksFromHub) = 5; % Local

    labelNames = {sprintf('hub (%u)',sum(sum(nodeLabels==1))),...
                  sprintf('feed-in (%u)',sum(sum(nodeLabels==2))),...
                  sprintf('feed-out (%u)',sum(sum(nodeLabels==3))),...
                  sprintf('feed-in-out (%u)',sum(sum(nodeLabels==4))),...
                  sprintf('local (%u)',sum(sum(nodeLabels==5)))};
end

end