% Originally TSQ_localcorrm, used in fetal heart rate analysis work.
% Takes in a distance matrix, and clusters it down to a reduced set
% Clusters down a set of operations based on their behavior
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-25
% Ben Fulcher, 2011-03-17
% ------------------------------------------------------------------------------

function [DistMat_cl,Cluster_Groupi,ord] = ClusterDown(DistMat,numClusters,opts)


if nargin < 2 || isempty(numClusters)
    numClusters = 5; % removes very highly-correlated operations 
end

% Additional parameters stored in the opts structure:
if nargin < 3
    opts = [];
end
if isfield(opts,'WhatDistance')
    % Input is a distance matrix based on correlation coefficient
    WhatDistance = opts.WhatDistance;
else
    WhatDistance = 'general';
end
if isfield(opts,'errth')
    errth = opts.errth;
else
    errth = 30; % set error threshold instead of just top number
end
if isfield(opts,'teststat') % necessary to do the thresholding
    teststat = opts.teststat;
end
if isfield(opts,'miorcorr')
    miorcorr = opts.miorcorr;
else
    miorcorr = 'corr'; % 'mi', 'corr'
end
if isfield(opts,'ClusterMeth')
    ClusterMeth = opts.ClusterMeth;
else
    ClusterMeth = 'linkage'; % 'linkage', 'kmedoids'
end
if isfield(opts,'LinkageMethod')
    LinkageMethod = opts.LinkageMethod;
else
    LinkageMethod = 'average';
end
if isfield(opts,'plotbig')
    plotbig = opts.plotbig;
else
    plotbig = 0;
end

% ------------------------------------------------------------------------------
% Cluster down to a reduced number of groups (some very highly-correlated operations):
% ------------------------------------------------------------------------------

figure('color','w');

% Do the linkage clustering:
fprintf(1,'Computing linkage information for %ux%u gene expression data...',length(DistMat),length(DistMat));
links = linkage(DistMat,LinkageMethod);
fprintf(1,' Done.\n');

% Get the dendrogram reordering:
subplot(5,1,1)
[h_dend,~,ord] = dendrogram(links,0);
% Reorder the distance matrix by dendrogram ordering: [could add optimalleaforder]
DistMat_cl = DistMat(ord,ord);   

% Make the dendrogram look aight
axis_pos = get(gca,'Position');
axis_pos(1) = 0.26;
axis_pos(3) = 0.42;
set(gca,'Position',axis_pos)
set(h_dend,'color','k','LineWidth',1)
set(gca,'XTickLabel',{})

% Compute clustering into groups for many different cluster numbers:
numClusterings = length(numClusters);
Cluster_Groupi = cell(numClusterings,1);
Cluster_Groupi_cl = cell(numClusterings,1);
ClusterCenters_cl = cell(numClusterings,1);
ClusterCenters = cell(numClusterings,1);

for i = 1:numClusterings
    fprintf(1,'Distance-based clustering with %u clusters\n',numClusters(i))
    
    % Cluster the dendrogram:
    T = cluster(links,'maxclust',numClusters(i));
    numClusters(i) = max(T);
    
    % Reorder members of each cluster by distance to other members of the cluster:
    Cluster_Groupi{i} = cell(numClusters(i),1);
    for j = 1:numClusters(i)
        Cluster_Groupi{i}{j} = find(T==j);
        if length(T==j) > 1
            % Sort by increasing sum of distances to other members of the cluster
            [~,ix] = sort(sum(DistMat(Cluster_Groupi{i}{j},Cluster_Groupi{i}{j})),'ascend');
            Cluster_Groupi{i}{j} = Cluster_Groupi{i}{j}(ix);
        end
    end

    % Select the closest to cluster centre in each group
    try
        ClusterCenters{i} = cellfun(@(x)x(1),Cluster_Groupi{i});

        Cluster_Groupi_cl{i} = cellfun(@(x) arrayfun(@(y)find(ord==y),x),Cluster_Groupi{i},'UniformOutput',0);
        ClusterCenters_cl{i} = arrayfun(@(y)find(ord==y),ClusterCenters{i});
    catch
        keyboard
    end
end


% ------------------------------------------------------------------------------
% Now plot it:
% ------------------------------------------------------------------------------
% Pick a given clustering and plot it
cLevel = min(1,numClusterings); % plot the first clustering

nlook = length(DistMat_cl);

% Plot as a similarity matrix:
subplot(5,1,2:5)
switch WhatDistance
case 'corr'
    % Input is a absolute correlation matrix:
    PlotCorrMat(1-DistMat_cl,1);
case 'general'
    % Input is a general distance matrix:
    PlotCorrMat(DistMat_cl,0);
end

% Add rectangles to indicate highly correlated clusters of statistics:
RectangleColors = BF_getcmap('accent',5,1);
for i = 1:numClusters
    % Cluster borders:
    rectangle('Position',[min(Cluster_Groupi_cl{cLevel}{i})-0.5,min(Cluster_Groupi_cl{cLevel}{i})-0.5, ...
                            length(Cluster_Groupi_cl{cLevel}{i}),length(Cluster_Groupi_cl{cLevel}{i})], ...
                    'EdgeColor',RectangleColors{1},'LineWidth',3);
                    
    % Cluster centers:
    rectangle('Position',[ClusterCenters_cl{cLevel}(i)-0.5,ClusterCenters_cl{cLevel}(i)-0.5,1,1], ...
                                'EdgeColor',RectangleColors{5},'LineWidth',3);
end


end