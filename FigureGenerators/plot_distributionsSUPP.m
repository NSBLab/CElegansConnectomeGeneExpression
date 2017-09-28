%===============================================================================
% SUPP distributions
%===============================================================================
% get data for reciprocally -unidirectionally and unconnected pairs
f = figure('color', 'w');
f.Position = [1000,200,1400,500];
[P1, S1, dataCell1] = coexpELCHRecUnidirUncon(C, G);
[dataCell2, P2, S2] = giveDistributions(C, G);

% combine data
dataCell = [dataCell1; dataCell2];
rgb_colorMatrix = GiveMeColors('AllDistributions');
colors = num2cell(rgb_colorMatrix, 2);
% extraParams = struct('customSpot','');
% extraParams = struct('theColors',colors);
extraParams.theColors = colors;
extraParams.customSpot = '.';
nEL = nansum(logical(dataCell{1}));
nCrec = nansum(logical(dataCell{2}));
nCuni = nansum(logical(dataCell{3}));
nUnc = nansum(logical(dataCell{4}));
nR = nansum(logical(dataCell{8}));
nFin = nansum(logical(dataCell{6}));
nFout = nansum(logical(dataCell{7}));
nP = nansum(logical(dataCell{5}));

JitteredParallelScatter(dataCell, true, true, false, extraParams)
ylabel('Gene coexpression, r_\phi', 'FontSize', 18);
set(gca,'Xtick', [1 2 3 4 5 6 7 8], 'XTickLabel',{sprintf('Electrical (%d pairs)', nEL), ...
                    sprintf('Chemical reciprocal (%d pairs)',nCrec),...
                    sprintf('Chemical unidirectional (%d pairs)',nCuni),...
                    sprintf('Unconnected (%d pairs)',nUnc),...
                    sprintf('Peripheral (%d pairs)',nP),...
                    sprintf('Feed-in (%d pairs)',nFin),...
                    sprintf('Feed-out (%d pairs)',nFout),...
                    sprintf('Rich (%d pairs)',nR)},...
                    'FontSize', 16);
xtickangle(20);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1]); %'Reciprocal', 'Unidirectional', 'Unconnected', 'Rich', 'Feed-in', 'Feed-out', 'Peripheral'}, 'FontSize', 16);
