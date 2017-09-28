function plotCoexpression(G,C,whatCorr,makeFigure,whatColors,sortNodesHow,kThreshold)

if nargin < 3
    whatCorr='';
end
if nargin < 4
    makeFigure = true;
end
if nargin < 5
    whatColors = 'ugly';
end
if nargin < 6
    sortNodesHow = 'xPos';
end
if nargin < 7
    D = GiveMeDefault();
    kThreshold = D.kHub;
end
%-------------------------------------------------------------------------------
% Sort neurons by degree (descending)
% Binary chemical+EM synapse network connectivity data
numNodes = length(GiveMeAdj(C));
switch sortNodesHow
case 'degree'
    [~,~, deg] = degrees_dir(GiveMeAdj(C));
    [sdeg,dind] = sort(deg,'descend');
    fprintf(1,'Gene coexpression ordered by degree\n');
case 'xPos'
    [~,dind] = sort(C.Pos(:,1),'descend'); % head-to-tail
    fprintf(1,'Gene coexpression ordered by x position\n');
case 'none'
    dind = 1:numNodes;
    fprintf(1,'Gene coexpression not sorted\n');

end

if makeFigure
    f = figure('color','w');
end

coExp = GiveMeCoexpression(G,whatCorr);
imagesc(coExp(dind,dind));

switch whatColors
case 'pretty'
    myBlues = BF_getcmap('blues',5,0);
    colormap([myBlues(1,:);1,1,1;BF_getcmap('reds',6)])
case 'ugly'
    colormap(flipud(BF_getcmap('spectral',9,0)))
end
caxis([-0.1,0.6])
axis square

% Threshold lines
% stopHere = find(sdeg < kThreshold,1,'first');
% hold on
% plot(0.5+[0,numNodes],ones(2,1)*stopHere,'k','LineWidth',1.5);
% plot(ones(2,1)*stopHere,0.5+[0,numNodes],'k','LineWidth',1.5);

end
