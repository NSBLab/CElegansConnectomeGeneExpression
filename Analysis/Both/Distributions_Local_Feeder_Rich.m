% ------------------------------------------------------------------------------
% Ben Fulcher, 2015-02-13
% Modified by Taqi Ali 24-7-15
% ------------------------------------------------------------------------------
% Simple plots for idiots, thanks ben
% ------------------------------------------------------------------------------

jitterOrMeanStd = 'meanstd'; % 'jitter'; 'meanstd';
kRich = 44;
whatAdj = 'rawBinary'; % 'rawBinary','raw weighted','zero binary','zero weighted'
networkType = 'bd'; % 'bd','bu','wd','wu'
whatNorm = 'Pearson_noLR'; % 'S1_Direct','S1_Indirect','S1_Direct_mod','S1_Indirect_mod' the different types of gene data sets
% ------------------------------------------------------------------------------

myColors = [BF_getcmap('spectral',4,1);BF_getcmap('set2',4,1);{zeros(1,3);ones(1,3)*0.4;ones(1,3)*0.6;ones(1,3)*0.8}];
myColors = myColors([1,3,4,9,10:12]);

linkedAdj = GiveMeAdj(C,whatAdj);
allCorrs = C.LineageDistance;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

maskTypes = {'rich','feeder','local','connected','connectedBoth','unconnected','unconnectedBoth'};
% maskTypes = {'connected','connectedBoth','connectedEither','unconnected','unconnectedBoth','unconnectedEither'};
numGroups = length(maskTypes);
CorrsCell = cell(numGroups,1);

for k = 1:numGroups
    mask = GiveMeMask(maskTypes{k},networkType,linkedAdj,{'degree',kRich},1);
    Corrs.(maskTypes{k}) = allCorrs(logical(mask.special) & ~isnan(allCorrs));
    CorrsCell{k} = Corrs.(maskTypes{k});
end

extraParams = struct;
extraParams.theColors = myColors;
extraParams.customSpot = '';
JitteredParallelScatter(CorrsCell,1,1,1,extraParams);
ax = gca; fig = gcf;
ax.XTick = 1:numGroups;
ax.XTickLabel = maskTypes;

set(gcf,'Position',[1817,474,392,218]);

% ------------------------------------------------------------------------------
% Compute p-values
% ------------------------------------------------------------------------------
for k1 = 1:numGroups
    for k2 = k1+1:numGroups
        [h,p,stats] = ttest2(CorrsCell{k1},CorrsCell{k2},'VarType','unequal');
        fprintf(1,'%s (%u)/%s (%u): p = %.2g\n',maskTypes{k1},length(CorrsCell{k1}),maskTypes{k2},length(CorrsCell{k2}),p);
    end
end

title(sprintf('Type %s and network %s',whatType,networkType));
