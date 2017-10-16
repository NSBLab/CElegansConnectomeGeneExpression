%===============================================================================
% MAIN TEXT:
%===============================================================================
% Plot coexpression between hub neurons that are born together/or early and non-hub
% neurons that are born together/or early.
%-------------------------------------------------------------------------------
relationship = 'born early';
D = GiveMeDefault();
khub = D.kHub;
early = 1000;
connected = 1;

[coexpEarly,ph,statsh, numRich, numNONRich] = compareBirthTimes(C,G, relationship, khub, early, connected);
extraParams = struct('customSpot','.');

rgb_colorMatrix = GiveMeColors('RichNONrich');
colors = num2cell(rgb_colorMatrix, 2);
% theColors1 = GiveMeColors('RichNONrich');
% theColors{1,:} = theColors1(1,:);
% theColors{2,:} = theColors1(2,:);
%theColors{1,:} = theColors1(2,:);
%theColors{2,:} = [.02 .56 .2];
extraParams.theColors = colors;
JitteredParallelScatter(coexpEarly,1,1,true,extraParams);


xLabels = {sprintf('Rich (%d) pairs',numRich), sprintf('Non-rich (%d) pairs', numNONRich)};


set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
set(gca,'box','off');
ylabel('Correlated gene expression, r_\phi','FontSize', 15);
% pos=get(sp3,'Position');
% set(sp3,'Position',[pos(1)*1.05, pos(2)*0.5, pos(3)*1, pos(4)*1.2]); % [left bottom width height]
