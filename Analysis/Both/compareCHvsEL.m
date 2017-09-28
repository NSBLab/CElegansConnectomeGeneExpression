% compare coexpression between electrical and chemical synapses
AdjCh = GiveMeAdj(C,'zeroBinary','ch');
AdjEl = GiveMeAdj(C,'zeroBinary','el');

AdjCh = triu(logical(AdjCh+AdjCh')+0);
AdjEl = triu(logical(AdjEl+AdjEl')+0);
NanMatrix = tril(nan(279,279));

AdjCh = AdjCh + NanMatrix; nCH = nansum(AdjCh(:));
AdjEl = AdjEl + NanMatrix; nEL = nansum(AdjEl(:));

CH = GiveMeCoexpression(G,'',true).*AdjCh;
EL = GiveMeCoexpression(G,'',true).*AdjEl;

CH = nonzeros(CH(:)); CH=CH(~isnan(CH)); dataCell{1} = CH;
EL = nonzeros(EL(:)); EL=EL(~isnan(EL)); dataCell{2} = EL;

[p,~,stats] = ranksum(CH,EL);
JitteredParallelScatter(dataCell);
ylabel('Gene coexpression, r_\phi', 'FontSize', 18);
set(gca,'Xtick', [1 2], 'XTickLabel',{sprintf('Chemical (%d pairs)', nCH), ...
    sprintf('Electrical (%d pairs)',nEL)}, 'FontSize', 16);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1]);
