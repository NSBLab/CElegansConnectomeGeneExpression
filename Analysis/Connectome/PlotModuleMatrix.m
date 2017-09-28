function PlotModuleMatrix(C,whatSolution,makeFigure)

%-------------------------------------------------------------------------------
if nargin < 2
    whatSolution = 'best';
end
if nargin < 3
    makeFigure = false;
end
%-------------------------------------------------------------------------------

gamma = 1;
tau = 0.1;
numRuns = 1000;
D = GiveMeDefault();
kHub = D.kHub;

%===============================================================================
% Load in data:
%===============================================================================
fileName = sprintf('ModuleAssignment_gamma%.2f_tau%.2f_runs%u.mat',gamma,tau,numRuns);
switch whatSolution
case 'ERMM'
    moduleAssignment = ERMMFit(C);
    load(fileName,'simMatrix');
    fprintf(1,'Getting module assignment from ERMM\n');
case {'best','consensus'}
    load(fileName)
    fprintf(1,'Loaded data from %s\n',fileName);
    % Get module assignments
    switch whatSolution
    case 'best'
        [~,isBest] = max(Q);
        moduleAssignment = Modules(:,isBest);
        fprintf(1,'Using optimal Q solution\n');
    case 'consensus'
        moduleAssignment = concModules;
        fprintf(1,'Using consensus solution\n');
    end
end
A = PlotMatrixRichFeederPeripheral(C,[],'',false,false);
%===============================================================================
ix = reorder_mod(simMatrix,moduleAssignment);
modAssOrd = moduleAssignment(ix);
A_ord = A(ix,ix);
[~,~,deg] = degrees_dir(simMatrix);
deg_ord = deg(ix);

% Position information:
posOrd = BF_NormalizeMatrix(C.Pos(ix,1),'maxmin');
anatomyType = LabelNeuronType(C,'anatomy');
anatomyTypeOrd = anatomyType(ix);
%-------------------------------------------------------------------------------

if makeFigure
    f = figure('color','w');
end
ax = gca; hold on
imagesc(A_ord);
axis('equal');
colors = GiveMeColors('richFeedInOutPeripheral');
colormap(colors);
set(gca,'Xtick',[],'XTickLabel',[]);
set(gca,'Ytick',[],'YTickLabel',[]);

%-------------------------------------------------------------------------------
% Add squares:
%-------------------------------------------------------------------------------
numModules = max(modAssOrd);
for i = 1:numModules
    moduleHere = [find(modAssOrd==i,1,'first'),find(modAssOrd==i,1,'last')];
    rectangle('Position',[moduleHere(1)-0.5,moduleHere(1)-0.5,diff(moduleHere)+1,diff(moduleHere)+1],'EdgeColor','w');
end

%-------------------------------------------------------------------------------
% Add rectangles:
%-------------------------------------------------------------------------------
neuronLabel = LabelNeuronType(C,'type');
neuronLabelOrd = neuronLabel(ix);
numNeurons = length(moduleAssignment);
rectThickness = 15;
colorsMod = BF_getcmap('set1',numModules,1);
[colorsNeuron,typeLabels] = GiveMeColors('InterneuronMotorSensoryMulti');
[colorsAnatomy,anatomyLabels] = GiveMeColors('anatomyType');
for j = 1:numNeurons
    % Module assignment
    rectangle('Position',[0.5-rectThickness,j-0.5,rectThickness,1], ...
                'FaceColor',colorsMod{modAssOrd(j)},'EdgeColor',colorsMod{modAssOrd(j)});

    % Neuron type
    labelj = (typeLabels==neuronLabelOrd(j));
    rectangle('Position',[0.5-rectThickness*2,j-0.5,rectThickness,1], ...
                'FaceColor',colorsNeuron(labelj,:),'EdgeColor',colorsNeuron(labelj,:));
    % Hub
    if deg_ord(j) > kHub
        c = 'r';
    else
        c = 'w';
    end
    rectangle('Position',[0.5-rectThickness*3,j-0.5,rectThickness,1], ...
                'FaceColor',c,'EdgeColor',c);

    % NeurotransmitterType
    % rectangle('Position',[0.5-rectThickness*3,j-0.5,rectThickness,1], ...
                % 'FaceColor',colors{modAssOrd(j)},'EdgeColor',colors{modAssOrd(j)});

    % Position:
    % rectangle('Position',[0.5-rectThickness*4,j-0.5,rectThickness,1], ...
    %             'FaceColor',ones(1,3)*posOrd(j),'EdgeColor',ones(1,3)*posOrd(j));

    % Head/tail/body:
    whatInd = (anatomyLabels==anatomyTypeOrd(j));
    rectangle('Position',[0.5-rectThickness*4,j-0.5,rectThickness,1], ...
                'FaceColor',colorsAnatomy(whatInd,:),'EdgeColor',colorsAnatomy(whatInd,:));

    % Boundaries:
    if j < numNeurons & modAssOrd(j+1)~=modAssOrd(j)
        plot([0.5-4*rectThickness,0.5],ones(2,1)*(j+0.5),'k','LineWidth',2)
    end
end

% Separators:
for k = 1:5
    plot(ones(2,1)*0.5-k*rectThickness,[0.5,numNeurons+0.5],'k')
end
plot([0.5-4*rectThickness,0.5],ones(2,1)*0.5,'k')
plot([0.5-4*rectThickness,0.5],ones(2,1)*(numNeurons+0.5),'k')

ax.XLim = [0.5 - 4*rectThickness,numNeurons+0.5];
ax.YLim = [0.5,numNeurons+0.5];
ax.XTick = 0.5+(-4*rectThickness+rectThickness/2:rectThickness:-rectThickness/2);
ax.XTickLabel = {'anatomyClass','Hub','NeuronType','Module'};
ax.XTickLabelRotation = 90;

end
