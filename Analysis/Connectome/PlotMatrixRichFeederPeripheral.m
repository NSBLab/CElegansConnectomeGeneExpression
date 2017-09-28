function A = PlotMatrixRichFeederPeripheral(C,kHub,edgeTypes,makeFigure,sortHow,addMeta)
%-------------------------------------------------------------------------------
% Plot a binary connectome, coloring edges by rich/feeder/peripheral
%-------------------------------------------------------------------------------

if nargin < 1
    C = load('CElegans_ConnectivityData.mat','C');
end
if nargin < 2 || isempty(kHub)
    D = GiveMeDefault();
    kHub = D.kHub;
    fprintf(1,'Hubs at k > 44, default\n');
end
if nargin < 3
    edgeTypes = '';
end
if nargin < 4
    makeFigure = true;
end
if nargin < 5
    sortHow = 'xPos';
end
if nargin < 6
    addMeta = true;
end

%-------------------------------------------------------------------------------

% Get adjacency matrix:
Adj = GiveMeAdj(C,'zeroBinary',edgeTypes);
[~,~,deg] = degrees_dir(Adj);
isHub = (deg > kHub);
fprintf(1,'%u hubs at k > %u\n',sum(isHub),kHub);

% Construct mask for edge types:
mask = nan(C.numNeurons,C.numNeurons);
mask(isHub,isHub) = 1; % rich
mask(isHub,~isHub) = 3; % feeder
mask(~isHub,isHub) = 2; % feeder
mask(~isHub,~isHub) = 4; % peripheral

% Generate masked adjacency matrix with labeled edge types, as A
A = logical(Adj).*mask;

% Sort nodes by degree?
switch sortHow
case 'degree'
    [~,ix] = sort(deg,'descend');
    A = A(ix,ix);
    fprintf(1,'Nodes sorted by degree (decreasing)\n');
case 'xPos'
    [~,ix] = sort(C.Pos(:,1),'ascend'); % head-to-tail
    A = A(ix,flipud(ix));
    fprintf(1,'Nodes sorted by x-position\n');
case 'none'
    ix = 1:length(Adj);
    fprintf(1,'Adj matrix nodes not sorted\n');
end

%-------------------------------------------------------------------------------
if makeFigure
    f = figure('color','w');
    f.Position = [1000         633         911         705];
end
ax = gca;
hold on
imagesc(A);
axis('equal');
colors = GiveMeColors('richFeedInOutPeripheral');
colormap(colors);
ylabel('Source neuron')
xlabel('Target neuron')
ax.XTick = [];
ax.YTick = [];

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Add rectangles:
%-------------------------------------------------------------------------------
if addMeta
    rectThickness = 25;
    numNeurons = length(A);

    neuronLabel = LabelNeuronType(C,'type');
    neuronLabelOrd = neuronLabel(ix);
    [colorsNeuron,typeLabels] = GiveMeColors('InterneuronMotorSensoryMulti');

    [~,ntTypes] = NeurotransmitterAnal(C,false);
    ntTypesOrd = ntTypes(ix);
    ntLabels = categories(ntTypes);
    colorsNeuroTrans = BF_getcmap('spectral',length(ntLabels),true);

    anatomyType = LabelNeuronType(C,'anatomy');
    anatomyTypeOrd = anatomyType(ix);
    [colorsAnatomy,anatomyLabels] = GiveMeColors('anatomyType');
    % posOrd = BF_NormalizeMatrix(C.Pos(ix,1),'maxmin');

    birthTimes = GiveMeBirthTimes(C);
    birthTimesOrd = BF_NormalizeMatrix(birthTimes(ix),'maxmin');

    degOrd = deg(ix);

    for j = 1:numNeurons
        % Neuron type
        labelj = (typeLabels==neuronLabelOrd(j));
        rectangle('Position',[0.5-rectThickness*1,j-0.5,rectThickness,1], ...
                    'FaceColor',colorsNeuron(labelj,:),'EdgeColor',colorsNeuron(labelj,:));
        % Hub
        if degOrd(j) > kHub
            c = 'r';
        else
            c = 'w';
        end
        rectangle('Position',[0.5-rectThickness*2,j-0.5,rectThickness,1], ...
                    'FaceColor',c,'EdgeColor',c);

        % NeurotransmitterType
        whatInd = (ntLabels==ntTypesOrd(j));
        rectangle('Position',[0.5-rectThickness*3,j-0.5,rectThickness,1], ...
                    'FaceColor',colorsNeuroTrans{whatInd},'EdgeColor',colorsNeuroTrans{whatInd});

        % Head/tail/body:
        whatInd = (anatomyLabels==anatomyTypeOrd(j));
        rectangle('Position',[0.5-rectThickness*4,j-0.5,rectThickness,1], ...
                    'FaceColor',colorsAnatomy(whatInd,:),'EdgeColor',colorsAnatomy(whatInd,:));

        % Birth time (earlier will be darker)
        rectangle('Position',[0.5-rectThickness*5,j-0.5,rectThickness,1], ...
                    'FaceColor',ones(1,3)*birthTimesOrd(j),'EdgeColor',ones(1,3)*birthTimesOrd(j));
    end

    % Separators:
    for k = 1:5
        plot(ones(2,1)*0.5-k*rectThickness,[0.5,numNeurons+0.5],'k')
    end
    plot([0.5-5*rectThickness,0.5],ones(2,1)*0.5,'k')
    plot([0.5-5*rectThickness,0.5],ones(2,1)*(numNeurons+0.5),'k')

    ax.XLim = [0.5 - 5*rectThickness,numNeurons+0.5];
    ax.YLim = [0.5,numNeurons+0.5];
    ax.XTick = 0.5+(-5*rectThickness+rectThickness/2:rectThickness:-rectThickness/2);
    ax.XTickLabel = {'Birth time','Anatomy','Neurotransmitter','Hub','Type'};
    ax.XTickLabelRotation = 90;
end


end
