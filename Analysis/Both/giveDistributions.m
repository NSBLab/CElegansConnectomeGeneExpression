% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute,
% make violin plots for rich/feeder/peripheral/unconnected for a chosen hub
% threshold. Can change to lineage matrix instead of coexpression matrix
% ------------------------------------------------------------------------------
function [dataCell, P2, S2] = giveDistributions(C, G, type)
    if nargin <3
        type = 'all';
    end
Adj = GiveMeAdj(C, 'zeroBinary', type);
[~,~,deg] = degrees_dir(Adj);
Adj = Adj|Adj';
Adj = triu(Adj,1);
NaNmatrix = nan(279,279);
NaNmatrix = tril(NaNmatrix,-1);
Adj = Adj+NaNmatrix;
correlation = GiveMeCoexpression(G,'',true); % C.LineageDistance;GiveMeCoexpression(G,'',true); %
if strcmp(type, 'ch')
    D = GiveMeDefault('synaptic');
elseif strcmp(type, 'all')
    D = GiveMeDefault('all');
end

k=D.kHub;

%[~,~,deg] = degrees_dir(Adj);
numNeurons = C.numNeurons;
hub = deg>k;
mask = zeros(numNeurons,numNeurons);
mask(hub, hub)=1;
mask(~hub, hub)=2;
mask(hub, ~hub)=3;
mask(~hub,~hub)=4;
mask = mask.*Adj;

rich = mask==1;
rich = nonzeros(rich.*correlation);
dataCell{4,1} = rich(~isnan(rich));

feedin = mask==2;
feedin = nonzeros(feedin.*correlation);
dataCell{2,1} = feedin(~isnan(feedin));


feedout = mask==3;
feedout = nonzeros(feedout.*correlation);
dataCell{3,1} = feedout(~isnan(feedout));

peripheral = mask==4;
peripheral = nonzeros(peripheral.*correlation);
dataCell{1,1} = peripheral(~isnan(peripheral));

%JitteredParallelScatterAA(dataCell)


[P2.RF,~,S2.RF] = ranksum(rich,feedin);
[P2.FinFout,~,S2.FinFout]  = ranksum(feedin, feedout); %,'Vartype','unequal','Tail','right');
[P2.FoutP,~,S2.FoutP]  = ranksum(feedout,peripheral);%,'Vartype','unequal','Tail','right');
[P2.FinP,~,S2.FinP]  = ranksum(feedin,peripheral); %,'Vartype','unequal','Tail','right');
[P2.RP,~,S2.RP]  = ranksum(rich,peripheral); %,'Vartype','unequal','Tail','right');
%[~,p6,~,stat6]  = ttest2(feedin,nonzeros(peripheral),'Vartype','unequal','Tail','right');
% [~,p1,~,stat1] = ttest2(rich,feedin,'Vartype','unequal','Tail','right');
% [~,p2,~,stat2]  = ttest2(feedin, feedout,'Vartype','unequal','Tail','right');
% [~,p3,~,stat3]  = ttest2(feedout,nonzeros(peripheral),'Vartype','unequal','Tail','right');
% [~,p4,~,stat4]  = ttest2(rich,feedout,'Vartype','unequal','Tail','right');
% [~,p5,~,stat5]  = ttest2(rich,nonzeros(peripheral),'Vartype','unequal','Tail','right');
% [~,p6,~,stat6]  = ttest2(feedin,nonzeros(peripheral),'Vartype','unequal','Tail','right');
end
