%===============================================================================
% MAIN TEXT:
%===============================================================================
% Plot coexpression as a function of degree for each neuron type.
% At each degree threshold we calculate degree coexpression between neurons where
% at least one neuron in a  pair is of a selected type and the degree of both neurons
% is higher than a threshold k.
%-------------------------------------------------------------------------------
% run calculations
%-------------------------------------------------------------------------------
[meanCoexp, numConnections, deg, types, stats] = typeDegree(C,G);

%===============================================================================
% Plot
%-------------------------------------------------------------------------------
f = figure('color', 'w');
f.Position = [1000,200,1000,500];
lineStyle = '-';
markerStyle = 'o';
% plot distributions for connections
numConnectionsLOG = log10(numConnections);
%meanCoexp(~any(~isnan(meanCoexp), 2),:)=[];
inds = find(isnan(meanCoexp(:,1)));
numConnectionsLOG(inds,:)= NaN;

colorList = GiveMeColors('InterneuronMotorSensoryMulti');
sp1 = subplot(2,5,1:3);
box off;
hold on
for t = 1:length(types)
    plot(1:max(deg),numConnectionsLOG(:,t), 'color', colorList(t,:),'LineWidth',2.5);
end
%H = bar(1:max(deg),log10(numConnections), 'stacked', 'EdgeColor', 'w'); xlim([1,60]);% ylim([0,3500]);
sp1.XTick = [];
pos = sp1.Position;
sp1.Position = [pos(1)*0.6, pos(2)*1, pos(3)*1, pos(4)*0.5]; % [left bottom width height]
ylabel('log_1_0(N_p_a_i_r_s)', 'FontSize', 15);%xlabel('strength>=s');
set(gca,'Ytick', [0 1 2 3 4], 'YTickLabel',[0 1 2 3 4], 'FontSize', 15);
%sp1.YTick = [0 1 2 3 4];
%sp1.YTickLabel = [0 1 2 3 4];
%sp1.FontSize = 12;
% for v=1:length(types)
%     set(H(v),'facecolor',colorList(v,:));
% end

%plot curves for each type
sp2 = subplot(2,5,6:8);
box off; hold on;
for i = 1:length(types)
    kr = 1:max(deg);
    plot(kr, meanCoexp(:,i), lineStyle, 'color',colorList(i,:),'LineWidth',2.5);
    %plot(1:max(deg), meanCoexp(:,i),markerStyle,'MarkerEdgeColor',colorList(i,:),...
    %'MarkerFaceColor',brighten(colorList(i,:),+0.5),'LineWidth',1,'MarkerSize',7); %xlim([0 60]);


pvalues = stats(:,i);
isSig = (pvalues < 0.05); % significantly higher than null

% mean (real data trajectory):
lineStyle = '-'; markerStyle = 'o';

if any(isSig)
    plot(kr(isSig),meanCoexp(isSig,i),markerStyle,'MarkerEdgeColor',colorList(i,:),...
        'MarkerFaceColor',brighten(colorList(i,:),+0.5),'LineWidth',1,'MarkerSize',6)
end
end

% mean trajectory

pos = sp2.Position;
sp2.Position = [pos(1)*0.6, pos(2)*2.5, pos(3)*1, pos(4)*0.8]; % [left bottom width height]
xlabel('Degree, k', 'FontSize', 15);
axisName = {'Median correlated', 'gene expression, r_\phi'};
ylabel(axisName, 'FontSize', 12);
set(gca,'Ytick', [0 .2 .4 .6], 'YTickLabel',[0 .2 .4 .6], 'FontSize', 15);
hold on;



sp3 = subplot(2,5,4:5);
[P,S,dataCell] = plotInterneuronsDegree(C,G);

rgb_colorMatrix = GiveMeColors('HubNonHub');
colors = num2cell(rgb_colorMatrix, 2);
extraParams = struct('customSpot','.');
extraParams.theColors = colors;
JitteredParallelScatter(dataCell,true,1,false, extraParams);
xLabels = {sprintf('Rich interneuron pairs (%d)',length(dataCell{1})), sprintf('Non-rich interneuron pairs (%d)', length(dataCell{2}))};
%xtickangle(20);
set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
set(gca,'box','off');
ylabel('Correlated gene expression, r_\phi','FontSize', 15);
hold on;

set(sp3, 'Position', [0.6 0.28 0.3 0.47]); 
%set(sp3,'Position',[pos(1), pos(2), pos(3), pos(4)]); % [left bottom width height]
%labels = {'A', 'B', 'C'}; 
%LabelCurrentAxes(labels); 