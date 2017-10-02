%===============================================================================
% Plot hub command interneurons and hub non-command interneurons
%===============================================================================

[dataCell, P,S, nC, nN] = plotInterneuronDistributions(C,G); 
JitteredParallelScatter(dataCell); 

xLabels = {sprintf('Hub command interneurons (%d)',nC), sprintf('Hub non-command interneurons (%d)', nN)};


set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
set(gca,'box','off');
ylabel('Gene coexpression, r_\phi','FontSize', 15);