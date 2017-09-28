function VisRichConns(C)
%-------------------------------------------------------------------------------
% Idea is to plot hub-hub connections through space
%-------------------------------------------------------------------------------

% labelWhat = 'hub';
labelWhat = 'type';
fontSize = 13;


% HEAD close-up
f = figure('color','w','Position',[318,529,1000,657]);
ax = cell(3,1);
ax{1} = subplot(2,4,1:2); hold on; box('on');
ax{2} = subplot(2,4,3:4); hold on; box('on');
ax{3} = subplot(2,4,5:8); hold on; box('on');
for i = 1:3
    ax{i}.FontSize = fontSize;
end
numAxes = length(ax);

% Neuron labeling
D = GiveMeDefault();
Adj = GiveMeAdj(C,'zeroBinary');
[~,~,deg] = degrees_dir(Adj);
isHub = (deg > D.kHub)';
switch labelWhat
case 'hub'
    neuronLabels = zeros(C.numNeurons,1);
    neuronLabels(isHub) = 1;
    neuronLabels = categorical(neuronLabels,[0,1],{'hub','nonhub'});
    theLabels = categories(neuronLabels);
    theColors = GiveMeColors('HubNonHub');
case 'type'
    neuronLabels = LabelNeuronType(C,'type');
    theLabels = {'interneuron','motor','sensory','multi'};
    theColors = GiveMeColors('InterneuronMotorSensoryMulti');
end
numLabels = length(theLabels);

xyPos = -C.Pos;

for i = 1:numAxes
    % Set axis to plot in
    axes(ax{i});

    %-------------------------------------------------------------------------------
    % Plot hub-hub connections
    %-------------------------------------------------------------------------------
    isRich = double(isHub)*double(isHub)';

    % Plot non-rich connections
    if i < 3
        isNotRichConn = ~isRich.*Adj;
        [theConns_i,theConns_j] = find(isNotRichConn);
        numConns = length(theConns_i);
        for j = 1:numConns
            plot([xyPos(theConns_i(j),1),xyPos(theConns_j(j),1)],[xyPos(theConns_i(j),2),xyPos(theConns_j(j),2)],...
                        '-','LineWidth',0.2,'color',ones(1,3)*0.85)
        end
    end

    % Plot rich connections
    isRichConn = isRich.*Adj;
    [theConns_i,theConns_j] = find(isRichConn);
    numConns = length(theConns_i);
    for j = 1:numConns
        plot([xyPos(theConns_i(j),1),xyPos(theConns_j(j),1)],[xyPos(theConns_i(j),2),xyPos(theConns_j(j),2)],...
                    'color',[1 .1 .07])
    end

    %-------------------------------------------------------------------------------
    % Plot neurons by hub/nonhub type
    %-------------------------------------------------------------------------------
    for k = 1:numLabels
        isType = neuronLabels==theLabels{k};
        for h = 1:2
            if h==1,
                r = isType & isHub;
                markerSize = 10;
                edgy = 'k';
            else
                r = isType & ~isHub;
                markerSize = 7;
                edgy = brighten(theColors(k,:),-0.4);
            end
            if i==3, markerSize = markerSize * 0.8; end
            plot(xyPos(r,1),xyPos(r,2),'o','MarkerSize',markerSize,...
                    'MarkerFaceColor',theColors(k,:),'MarkerEdgeColor',edgy);
        end
    end
    %-------------------------------------------------------------------------------

    %-------------------------------------------------------------------------------
    % Now adjust axis limits:
    %-------------------------------------------------------------------------------
    axis('equal')
    switch i
    case 1 % tail
        ax{i}.XLim = sort(-[-1.03,-0.93]);
        ax{i}.YLim = sort(-[-0.053,-0.013]);
        title('tail')
        ax{i}.Position = [0.5422    0.3388    0.3628    0.3412];
    case 2 % head
        ax{i}.XLim = sort(-[0,0.142]);
        ax{i}.YLim = sort(-[-0.025,0.03]);
        ax{i}.Position = [0.1300    0.3388    0.3628    0.3412];
        title('head')
    case 3 % all
        for k = 1:2
            rectangle('Position',[ax{k}.XLim(1),ax{k}.YLim(1),diff(ax{k}.XLim),diff(ax{k}.YLim)],...
                        'EdgeColor','k','LineStyle',':','LineWidth',2)
        end
        ax{i}.XLim = [min(xyPos(:,1))-0.02,max(xyPos(:,1))+0.02];
        ax{i}.YLim = [min(xyPos(:,2))-0.02,max(xyPos(:,2))+0.02];
        ax{i}.Position = [0.1300    0.0963    0.7750    0.3412];
    end
    xlabel(ax{i},'anterior-posterior axis')
    ylabel(ax{i},'dorsal-ventral axis')
end

end
