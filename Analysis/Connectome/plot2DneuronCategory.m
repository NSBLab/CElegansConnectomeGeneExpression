%-------------------------------------------------------------------------------
% Aurina Arnatkeviciute
% plot hierarchical annotations (head/tail/other) in 2D space
% colorcode according to hierarchical annotation

%-------------------------------------------------------------------------------%%
% head and tail
%-------------------------------------------------------------------------------
position = C.Pos;

figure('color','w');
ax1 = subplot(3,1,1); hold on

anatomyTypes = LabelNeuronType(C,'anatomy');
anatomyTypeLabels = {'head','body','tail'};
numAnatomyTypes = 3;
[c,clabels] = GiveMeColors('directedAnatomy');
c = c(cellfun(@(x)find(strcmp(clabels,x)),{'HeadHead','BodyBody','TailTail'}),:);
for j = 1:numAnatomyTypes
    isType = anatomyTypes==anatomyTypeLabels{j};
    plot(position(isType,1),position(isType,2),'o','MarkerSize',7,...
                'MarkerFaceColor',c(j,:),'MarkerEdgeColor',brighten(c(j,:),-0.4));
end
legend(anatomyTypeLabels)
axis('equal')
ax1.XLim = [min(position(:,1)),max(position(:,1))];
ax1.YLim = [min(position(:,2)),max(position(:,2))];

%-------------------------------------------------------------------------------
%% sensory/motor/interneuron
%-------------------------------------------------------------------------------
systemTypes = LabelNeuronType(C);
systemTypeLabels = {'interneuron','motor','sensory','multi'};
numSystemTypes = length(systemTypeLabels);
c = GiveMeColors('InterneuronMotorSensoryMulti');

ax2 = subplot(3,1,2); hold on
for j = 1:numSystemTypes
    isType = systemTypes==systemTypeLabels{j};
    plot(position(isType,1),position(isType,2),'o','MarkerSize',7,...
                'MarkerFaceColor',c(j,:),'MarkerEdgeColor',brighten(c(j,:),-0.4));
end
legend(systemTypeLabels)
axis('equal')
ax2.XLim = [min(position(:,1)),max(position(:,1))];
ax2.YLim = [min(position(:,2)),max(position(:,2))];
% colormap(jet);
% colorbar; %caxis([0 25]);

%-------------------------------------------------------------------------------
%% Neurotransmitter types
%-------------------------------------------------------------------------------
[~,ntTypes] = NeurotransmitterAnal(C,false);
ntLabels = categories(ntTypes);
numNTLabels = length(ntLabels);
c = BF_getcmap('spectral',numNTLabels,0);
ax3 = subplot(3,1,3); hold on
for j = 1:numNTLabels
    isType = ntTypes==ntLabels{j};
    plot(position(isType,1),position(isType,2),'o','MarkerSize',7,...
                'MarkerFaceColor',c(j,:),'MarkerEdgeColor','k');
end
legend(ntLabels)
axis('equal')
ax3.XLim = [min(position(:,1)),max(position(:,1))];
ax3.YLim = [min(position(:,2)),max(position(:,2))];
