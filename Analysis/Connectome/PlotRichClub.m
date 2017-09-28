function PlotRichClub(synapseType,networkType,nullModel,numIter,numNullNetworks)
%-------------------------------------------------------------------------------
% Loads in data from a computation using networkRC,
% Plots RC curve and degree distribution
%-------------------------------------------------------------------------------

doRankSum = false; % Welch's t-test or ranksum test?
meanOrMedian = 'mean'; % mean distance

%-------------------------------------------------------------------------------
% Load in the data
%-------------------------------------------------------------------------------
fileNameIn = fullfile('Data',sprintf('richClubResults_%s-%s-%s_%u-%uiter.mat',synapseType,...
                                            networkType,nullModel, numNullNetworks,numIter));
load(fileNameIn);
fprintf(1,'Loaded %s...\n',fileNameIn);
%-------------------------------------------------------------------------------

f = figure('color','w');
axMain = subplot(3,1,2:3); hold on;

%===============================================================================
% Distance curve
%===============================================================================
distData = GiveMeDist(C);
kRange = min(deg):max(deg);
% Add a baseline (1)
plot([min(deg),max(deg)],ones(2,1),':k');
dMean = zeros(length(kRange),2);
dAll = cell(length(kRange),2);
pVals = zeros(length(kRange),1);
for i = 1:length(kRange)
    isHub = double(deg > kRange(i));
    isRich = isHub'*isHub;
    isRich(eye(logical(size(isRich)))) = 0;
    notRich = ~isRich;
    notRich(eye(logical(size(isRich)))) = 0;
    dAll{i,1} = distData(Adj & isRich);
    dAll{i,2} = distData(Adj & notRich);

    % Difference
    if all(isnan(dAll{i,1})) || all(isnan(dAll{i,2}))
        pVals(i) = NaN;
    else
        if doRankSum
            pVals(i) = ranksum(dAll{i,1},dAll{i,2},'tail','right');
        else
            [~,pVals(i)] = ttest2(dAll{i,1},dAll{i,2},'Tail','right','Vartype','unequal');
        end
    end
end

switch meanOrMedian
case 'mean'
    dMiddle = cellfun(@mean,dAll);
case 'median'
    dMiddle = cellfun(@median,dAll);
end
dMiddleNorm = 1+0.4*(dMiddle(:,1)-dMiddle(1,1))/(max(dMiddle(:,1))-dMiddle(1,1));
h_dCurve = plot(kRange,dMiddleNorm,'color',[.58 .34 .92],'LineWidth',2);
sigD = pVals < 0.05;
plot(kRange(sigD),dMiddleNorm(sigD),'o','MarkerEdgeColor',[.58 .34 .92],...
        'MarkerFaceColor',brighten([.58 .34 .92],0.5),'LineWidth',1,'MarkerSize',5.3)

%===============================================================================
% Add rich-club curve:
%===============================================================================
h_RCcurve = plot(PhiNormMean, '-','Color','r','LineWidth',2);
% Significance at p = 0.05 as circles
plot(find(isSig),PhiNormMean(isSig),'o','MarkerEdgeColor','r',...
        'MarkerFaceColor',[1 .41 .38],'LineWidth',1,'MarkerSize',5.3)

% Add shaded rectangle
ylimNow = [axMain.YLim(1),axMain.YLim(2)];
h_rect = rectangle('Position',[43.5,ylimNow(1),20,ylimNow(2)-ylimNow(1)+.1],...
                    'EdgeColor','none','FaceColor',ones(3,1)*0.90);
uistack(h_rect,'bottom');
axMain.YLim = ylimNow;

% Set axis limits, legend, labels:
axMain.YLim = [0.9 max([max(PhiNormMean),max(dMiddleNorm)])+0.1];
xlabel('Degree, k','fontsize',14);
ylabel('\Phi_{norm}, d_{norm}','fontsize',14);
get(gca, 'XTick');
set(gca, 'FontSize', 12)
box off;
legend([h_RCcurve,h_dCurve],{'\Phi_{norm}','d_{norm}'},'Location','NorthWest','fontsize',12);
legend('boxoff') 
box('off')
axMain.XLim = [min(deg)-0.5,max(deg)+2];
LabelCurrentAxes('B',gca,26,'topRight')
%yyaxis right
%axMain.YLim = [0.9 max([max(PhiNormMean),max(dMiddleNorm)])+0.1];
%ylabel('\d_{Euclnorm}');

%===============================================================================
% Color degree distribution by different neuron properties
%===============================================================================
kUnique = unique(deg);
colormap(GiveMeColors('InterneuronMotorSensoryMulti'));

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
ylabel('Frequency','fontsize',14);
get(gca, 'YTick');
set(gca, 'FontSize', 12)
legend(theLabelsLegend,'Position',[0.413,0.781,0.176,0.072], 'fontsize',12); %'Location','NorthEast')
legend('boxoff') 
box('off')
axDD_type.Position = [0.1300    0.6384    0.7750    0.2823];
LabelCurrentAxes('A',gca,26,'topRight')

f.Position = [738,715,560,353];

end
