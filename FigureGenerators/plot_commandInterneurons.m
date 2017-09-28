% plot distributions for hub command interneurons (10) and hub neurons that are
% not command interneurons (6). 
% take pairs, not connections. 

[dataCellALL, P,S] = plotInterneuronDistributions(C,G); 


dataCell = dataCellALL(3:4); 

rgb_colorMatrix = GiveMeColors('HubNonHub');
colors = num2cell(rgb_colorMatrix, 2);
extraParams = struct('customSpot','.');
extraParams.theColors = colors;
JitteredParallelScatter(dataCell,true,1,true, extraParams);
xLabels = {'Command interneurons', 'Non-command interneurons'};
%xtickangle(20);
set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
set(gca,'box','off');
ylabel('Gene coexpression, r_\phi','FontSize', 15);