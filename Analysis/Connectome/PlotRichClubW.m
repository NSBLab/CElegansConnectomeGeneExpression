function h_RCcurveW = PlotRichClubW(synapseType,networkType,nullModel, numIter,numNullNetworks, whatColor, makefigure)
%-------------------------------------------------------------------------------
% Loads in data from a computation using networkRC,
% Plots RC curve and degree distribution
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Load in the data
%-------------------------------------------------------------------------------
fileNameIn = fullfile('Data',sprintf('richClubResults_%s-%s-%s_%u-%uiter.mat',synapseType,...
                                            networkType,nullModel, numNullNetworks,numIter));
load(fileNameIn);
fprintf(1,'Loaded %s...\n',fileNameIn);

%-------------------------------------------------------------------------------
if makefigure
f = figure('color','w');
end
axMain = subplot(3,1,2:3); hold on;

%===============================================================================
% Distance curve
%===============================================================================
% distData = GiveMeDist(C);
% kRange = min(deg):max(deg);
% % Add a baseline (1)
% plot([min(deg),max(deg)],ones(2,1),':k');
% dMean = zeros(length(kRange),2);
% dAll = cell(length(kRange),2);
% pVals = zeros(length(kRange),1);
% for i = 1:length(kRange)
%     isHub = double(deg > kRange(i));
%     isRich = isHub'*isHub;
%     isRich(eye(logical(size(isRich)))) = 0;
%     notRich = ~isRich;
%     notRich(eye(logical(size(isRich)))) = 0;
%     dAll{i,1} = distData(Adj & isRich);
%     dAll{i,2} = distData(Adj & notRich);
%
%     % Difference:
%     [~,pVals(i)] = ttest2(dAll{i,1},dAll{i,2},'Tail','right','Vartype','unequal'); %,'Tail','right');
% end
% dMean = cellfun(@mean,dAll);
% dMeanNorm = 1+0.8*(dMean(:,1)-dMean(1,1))/(max(dMean(:,1))-dMean(1,1));
% h_dCurve = plot(kRange,dMeanNorm,'color',[.96 .78 .1],'LineWidth',2);
% sigD = pVals < 0.05;
% plot(kRange(sigD),dMeanNorm(sigD),'o','MarkerEdgeColor',[.96 .78 .1],...
%         'MarkerFaceColor',brighten([.96 .78 .1],0.5),'LineWidth',1,'MarkerSize',5.3)

%===============================================================================
% Add rich-club curve:
%===============================================================================
h_RCcurveW = plot(PhiNormMean, '-','Color',whatColor(1,:),'LineWidth',2);
% Significance at p = 0.05 as circles
plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor',whatColor(1,:),...
        'MarkerFaceColor',whatColor(2,:),'LineWidth',1,'MarkerSize',5.3)

% Add shaded rectangle
ylimNow = [axMain.YLim(1),axMain.YLim(2)];
h_rect = rectangle('Position',[43.5,ylimNow(1),20,2.3],...
                    'EdgeColor','none','FaceColor',ones(3,1)*0.90);
uistack(h_rect,'bottom');
axMain.YLim = ylimNow;

% Set axis limits, legend, labels:
axMain.YLim = [0.9 max([max(PhiNormMean),max(PhiNormMean)])+0.2];
xlabel('Degree, k');
ylabel('\Phi_{norm}');
box off;

axMain.XLim = [min(deg)-0.5,max(deg)+2];
plot([min(deg),max(deg)],ones(2,1),':k');
%LabelCurrentAxes('B',gca,26,'topRight')

hold on;

%===============================================================================
% ---DEGREE DISTRIBUTION---
%===============================================================================
% axDD = subplot(5,1,2);
% ~~~VANILLA DEGREE DISTRIBUTION~~~
% histogram(deg, 100,'EdgeColor',[0 0 0],'FaceColor',[.45 .45 .45]);
% axDD.XLim = axMain.XLim;
% axDD.XTickLabel = {};
% box off;

%===============================================================================
% Color degree distribution by different neuron properties
%===============================================================================
% axDD_NT = subplot(5,1,2);
% [k,NTLabels] = NeurotransmitterAnal(C,false);
% theLabels = categories(NTLabels);
% numLabels = length(theLabels);
% % Make a categorical matrix matched to unique degrees:
kUnique = unique(deg);
% typeMatrix = zeros(length(kUnique),numLabels);
% for i = 1:length(kUnique)
%     for j = 1:numLabels
%         typeMatrix(i,j) = sum(deg(NTLabels==theLabels{j})==kUnique(i));
%     end
% end
% bar(kUnique,typeMatrix,'stacked')
colormap(GiveMeColors('InterneuronMotorSensoryMulti'));
% legend(theLabels,'Location','best','FontSize',8)
% axDD_NT.XLim = axMain.XLim;
% axDD_NT.XTickLabel = {};
% box off;

%===============================================================================
% Now do neuron types
%===============================================================================
neuronLabel = LabelNeuronType(C);

% Now make a typeMatrix (for neurotransmitters)
% theLabels = categories(neuronLabel);
theLabels = {'interneuron','motor','sensory','multi'};
theLabelsLegend = {sprintf('Interneuron (%u)',sum(neuronLabel=='interneuron')),...
            sprintf('Motor (%u)',sum(neuronLabel=='motor')),...
            sprintf('Sensory (%u)',sum(neuronLabel=='sensory')),...
            sprintf('Multiple (%u)',sum(neuronLabel=='multi'))};
numLabels = length(theLabels);
typeMatrix_NT = zeros(length(kUnique),numLabels);
for i = 1:length(kUnique)
    for j = 1:numLabels
        typeMatrix_NT(i,j) = sum(deg(neuronLabel==theLabels{j})==kUnique(i));
    end
end

axDD_type = subplot(3,1,1);
bar(kUnique,typeMatrix_NT,'stacked','LineWidth',0.2)
axDD_type.XLim = axMain.XLim;
axDD_type.XTickLabel = {};
ylabel('Frequency');
legend(theLabelsLegend,'Position',[0.413,0.781,0.176,0.072]); %'Location','NorthEast')
box('off')
axDD_type.Position = [0.1300    0.6384    0.7750    0.2823];
%LabelCurrentAxes('A',gca,26,'topRight')

f.Position = [738,715,560,353];

end
