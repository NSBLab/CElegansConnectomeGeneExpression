function [A,dist,ord] = plotGeneExpression(G,C,doInterneurons,sortGenesHow,extraHistogram)
%% reorders gene expression data according to clusters
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute
% ------------------------------------------------------------------------------
% Reorders gene expression data based on binary distance between vectors
% ------------------------------------------------------------------------------

if nargin < 3
    doInterneurons = false;
end
if nargin < 4
    sortGenesHow = 'descend';
end
if nargin < 5
    extraHistogram = false;
end
%-------------------------------------------------------------------------------

% Gene expression data
expData = (G.GeneExpData.Direct);
[numNeurons,numGenes] = size(expData);
[~,~, deg] = degrees_dir(GiveMeAdj(C));

% Only look at the subset of interneurons (but full connectivity data)
if doInterneurons
    fprintf(1,'Looking just at interneurons\n');
    isInterneuron = C.RegionM(:,10);
    expData = expData(isInterneuron,:);
    deg = deg(isInterneuron);
end

%-------------------------------------------------------------------------------
% Reorder genes
switch sortGenesHow
case 'descend'
    % Sort genes by number of expressed neurons:
    [~,ord] = sort(sum(expData),'descend');
case 'cluster'
    % Cluster by a gene expression similarity measure:
    % dist = 1-GiveMeCoexpression(G); % by binary pearsons
    dist = 'euclidean'; % by euclidean
    ord = BF_ClusterReorder(expData',dist);
end

%-------------------------------------------------------------------------------
% Sort neurons by degree (descending)
% Binary chemical synapse network connectivity data
% [sdeg,dind] = sort(deg,'descend');
[~,dind] = sort(C.Pos(:,1),'descend'); % head-to-tail
fprintf(1,'Gene expression ordered by x position\n');

% Reorder expression matrix
expData_ord = expData(dind,ord);

%-------------------------------------------------------------------------------
if extraHistogram
    f = figure('color','w');
end

% Set colormap:
colors = BF_getcmap('spectral',5,1);
col = [0 0 0; colors{4}];
colormap(col);

% Bar
if extraHistogram
    subplot(5,1,1); ax1 = gca;
    bar(sum(expData_ord,1),'EdgeColor',colors{4},'FaceColor','w')
    ax1.XLim = [0.5000,numGenes+0.5000];
    ax1.XTick = [];
    ax1.YLim = [0,max(sum(expData_ord,1))];
    ax1.Position = [0.1300    0.7548    0.7750    0.1243];
    ylabel('Number neurons')
end

% Colorplot
if extraHistogram
    subplot(5,1,2:5);
end
ax = gca;
imagesc(expData_ord);
xlabel('948 genes');
ylabel('279 neurons');
ax.YTick = [];
ax.XTick = [];
ax.YTickLabel = [];
ax.XTickLabel = [];
ax.XLim = [0.5000,numGenes+0.5000];
if extraHistogram
    f.Position = [997   664   609   349];
end

end
