%-------------------------------------------------------------------------------
% Do neurons that connect to hubs (or receive connections from hubs) display
% different properties to each other?
%-------------------------------------------------------------------------------

D = GiveMeDefault();

%-------------------------------------------------------------------------------
% Connectivity data:
Adj = GiveMeAdj(C,'zeroBinary');
[~,~,k] = degrees_dir(Adj);
numNeurons = length(Adj);

%-------------------------------------------------------------------------------
isHub = double(k > D.kHub);
isHubConn = logical((isHub'*isHub | isHub'*~isHub | ~isHub'*isHub) & Adj);
% isRichConn = Adj & isRich;

% First classify nodes
nodeTypes = struct();
nodeTypes.hub = logical(isHub)';
nodeTypes.feedIn = ~nodeTypes.hub & any(isHubConn,2); % any targets of a non-hub source node are hubs
nodeTypes.feedOut = (~nodeTypes.hub' & any(isHubConn,1))'; % any sources of a non-hub target node are hubs
nodeTypes.feedBoth = nodeTypes.feedIn & nodeTypes.feedOut; % project to and receive projections from hubs
nodeTypes.feedIn(nodeTypes.feedBoth) = false; % can't be both feed-in and feed-out
nodeTypes.feedOut(nodeTypes.feedBoth) = false; % can't be both feed-in and feed-out
nodeTypes.peripheral = ~nodeTypes.hub & ~nodeTypes.feedOut & ~nodeTypes.feedIn; % all other nodes are labeled 'peripheral'
nodeTypeLabels = fieldnames(nodeTypes);
numNodeTypes = length(nodeTypeLabels);


f = figure('color','w'); ax = gca;

%===============================================================================
% Neurotransmitter types:
%===============================================================================
[~,ntTypes] = NeurotransmitterAnal(C,false);
ntLabels = categories(ntTypes);
numNTLabels = length(ntLabels);

identities = zeros(numNodeTypes,numNTLabels);
for i = 1:numNodeTypes
    isNodeType = nodeTypes.(nodeTypeLabels{i});
    for j = 1:numNTLabels
        identities(i,j) = sum(ntTypes(isNodeType)==ntLabels{j});% /sum(isNodeType);
        if ismember(nodeTypeLabels{i},{'hub','feedIn','feedOut'})
            fprintf(1,'%s nodes are %.2f%% %s\n',nodeTypeLabels{i},...
                            100*identities(i,j)/sum(isNodeType),ntLabels{j});
        end
    end
end
f = figure('color','w');
ax = gca;
bar(identities,'stacked')
ax.XTick = 1:5;
ax.XTickLabel = nodeTypeLabels;
colormap(BF_getcmap('spectral',numNTLabels,0))
cB = colorbar;
cB.Ticks = linspace(1.5,8.5,numNTLabels);
cB.TickLabels = ntLabels;
ax.Position = [0.1184    0.1100    0.6637    0.8150];
ylabel('Number of neurons')
title('Neurotransmitter system of different neuron types')

%-------------------------------------------------------------------------------
% Get anatomy identities:
%-------------------------------------------------------------------------------
anatomyTypes = LabelNeuronType(C,'anatomy');
anatomyTypeLabels = {'head','body','tail'};
numAnatomyTypes = 3;

identities = zeros(numNodeTypes,numAnatomyTypes);
for i = 1:numNodeTypes
    isNodeType = nodeTypes.(nodeTypeLabels{i});
    for j = 1:numAnatomyTypes
        identities(i,j) = sum(anatomyTypes(isNodeType)==anatomyTypeLabels{j});% /sum(isNodeType);
        if ismember(nodeTypeLabels{i},{'feedIn','feedOut'})
            fprintf(1,'%s nodes are %.2f%% %s\n',nodeTypeLabels{i},...
                            100*identities(i,j)/sum(isNodeType),anatomyTypeLabels{j});
        end
    end
end
f = figure('color','w'); axAnatomy = gca;
bar(identities,'stacked')
axAnatomy.XTick = 1:5;
axAnatomy.XTickLabel = nodeTypeLabels;
[rgb_colorMatrix,labels] = GiveMeColors('directedAnatomy');
colormap(rgb_colorMatrix(cellfun(@(x)find(strcmp(labels,x)),{'HeadHead','BodyBody','TailTail'}),:))
cB = colorbar;
cB.Ticks = linspace(1+1/3,3-1/3,3);
cB.TickLabels = {'Head','Body','Tail'};
axAnatomy.Position = [0.1184    0.1100    0.6637    0.8150];
ylabel('Number of neurons')

%-------------------------------------------------------------------------------
% Get system type identities:
%-------------------------------------------------------------------------------
systemTypes = LabelNeuronType(C,'type');
systemTypeLabels = {'interneuron','motor','sensory','multi'};
numSystemTypes = length(systemTypeLabels);

identities = zeros(numNodeTypes,numSystemTypes);
for i = 1:numNodeTypes
    isNodeType = nodeTypes.(nodeTypeLabels{i});
    for j = 1:numSystemTypes
        identities(i,j) = sum(systemTypes(isNodeType)==systemTypeLabels{j});% /sum(isNodeType);
        if ismember(nodeTypeLabels{i},{'feedIn','feedOut'})
            fprintf(1,'%s nodes are %.2f%% %s\n',nodeTypeLabels{i},...
                            100*identities(i,j)/sum(isNodeType),systemTypeLabels{j});
        end
    end
end
f = figure('color','w'); axSystem = gca;
bar(identities,'stacked')
axSystem.XTick = 1:5;
axSystem.XTickLabel = nodeTypeLabels;
colormap(GiveMeColors('InterneuronMotorSensoryMulti'))
cB = colorbar;
cB.Ticks = linspace(1+1/3,4-1/3,4);
cB.TickLabels = systemTypeLabels;
axSystem.Position = [0.1184    0.1100    0.6637    0.8150];
ylabel('Number of neurons')
% SaveAllFigures('feedInOutNodes','eps',false);

%-------------------------------------------------------------------------------
% Birth times?
%-------------------------------------------------------------------------------
birthTimes = GiveMeBirthTimes(C);
dataCell = cell(5,1);
for i = 1:5
    dataCell{i} = birthTimes(nodeTypes.(nodeTypeLabels{i}));
end
JitteredParallelScatter(dataCell,true,true,true);
ax = gca;
ax.XTick = 1:5;
ax.XTickLabel = nodeTypeLabels;
ylabel('Birth time')

propEarly = zeros(5,1);
for i = 1:5
    propEarly(i) = mean(birthTimes(nodeTypes.(nodeTypeLabels{i})) < 1000);
end
f = figure('color','w'); ax = gca;
bar(propEarly)
ax.XTick = 1:5;
ax.XTickLabel = nodeTypeLabels;
