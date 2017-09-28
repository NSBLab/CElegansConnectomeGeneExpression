%% [p, stats] = coexpConnUncon (C, G, connectivityType, anatomyPart, dorankedTest)
% ------------------------------------------------------------------------------
% function compares and plots coexpression between pairs of neurons that are connected with electrical, chemical synapses as well as pairs of unconnected neurons.
% links in Celegans connectome.
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure), G (gene structure) should be loaded;
% coexpMeasure - choose coexpression measure as defined in G.Corr. default Pearson_noLR
% dorankedTest - true (will do ranksum test to compare distributions), false - will do ttest2 to compare distributions.
%% ------------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------------
% p - p value from the test comparing distributions for connected and unconnected links
% stats - stats output from the test comparing distributions for connected and unconnected links.
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-12-20
% ------------------------------------------------------------------------------
function [dataCell, S,P] = coexpELCHUncon (C, G, dorankedTest, doOnly)

if nargin < 3
    dorankedTest = true;
    fprintf(1,'Using rank test BY DEFAULT\n');
end

if nargin < 4 
    doOnly = false;
end

Adj = GiveMeAdj(C); 
maskisLR = GiveMeLRMask(C);

if doOnly == true
    mask = GiveMeMask(C,'hubType','bu',Adj);
    Chemical = mask.chemicalOnly & ~maskisLR;
    Electrical = mask.electricalOnly & ~maskisLR;

    %Unc = mask.none;
else
    mask = GiveMeMask(C,'hubType','bu', Adj);
    Chemical=GiveMeAdj(C,'zeroBinary', 'ch') & ~maskisLR;
    Electrical=GiveMeAdj(C,'zeroBinary', 'el') & ~maskisLR;
    
end

coexpression = GiveMeCoexpression(G,[],false);
coexpression = coexpression.*~maskisLR; 

Both = mask.chelboth|mask.chelboth' & ~maskisLR;

Con = mask.either;
Con = Con|Con' & ~maskisLR;

Unc = ~Con & ~maskisLR; 
Con = triu(Con,1);
Unc = triu(Unc,1); 


Uncon = coexpression.*Unc;
Chemical = Chemical|Chemical'; 
Chemical = triu(coexpression.*Chemical,1);
Electrical = Electrical|Electrical'; 
Electrical = triu(coexpression.*Electrical,1);
Con = coexpression.*Con;
Both = triu(coexpression.*Both,1);


% exclude zeros
Uncon = nonzeros(Uncon);

Chemical = nonzeros(Chemical);

Electrical = nonzeros(Electrical);

Both = nonzeros(Both);

Con = nonzeros(Con);


% compare distributions
if dorankedTest==false
    [~,P.pEC,~,S.statsEC] = ttest2(Electrical,Chemical, 'Vartype','unequal');
    [~,P.pCU,~,S.statsCU] = ttest2(Con,Uncon, 'Vartype','unequal');
elseif dorankedTest==true
    [P.pEC,~,S.statsEC] = ranksum(Electrical,Chemical);
    [P.pCU,~,S.statsCU] = ranksum(Con, Uncon);
    %[P.pBC,~,S.statsBC] = ranksum(Both, Chemical);
    %[P.pBE,~,S.statsBE] = ranksum(Both, Electrical);
    
end

if doOnly==true
    
    dataCell{1,1} = Electrical; dataCell{2,1} = Chemical; dataCell{3,1} = Both; dataCell{4,1} = Uncon;
else
    dataCell{1,1} = Electrical; dataCell{2,1} = Chemical; dataCell{3,1} = Uncon;
end


% xLabels = {'Connected', 'Unconnected'};
%
% JitteredParallelScatter(dataCell)
%
% set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% set(gca,'box','off');
% ylabel('Coexpression','FontSize', 15);
end
