%% [p, stats] = coexpConnUncon (C, G, connectivityType, anatomyPart, dorankedTest)
% ------------------------------------------------------------------------------
% function compares and plots coexpression for connected and unconnected
% links in Celegans connectome.
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure), G (gene structure) should be loaded;
% coexpMeasure - choose coexpression measure as defined in G.Corr. default Pearson_noLR
% connectivityType: chemical synapses - 'ch', electrical synapses - 'el', chemical + electrical synapses - 'all'
% anatomyPart: perform calculation on the whole nervous system - 'wholeNervousSystem'; perform calculations only on head neurons - 'headOnly'
% dorankedTest - true (will do ranksum test to compare distributions), false - will do ttest2 to compare distributions.
%% ------------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------------
% p - p value from the test comparing distributions for connected and unconnected links
% stats - stats output from the test comparing distributions for connected and unconnected links.
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-12-20
% ------------------------------------------------------------------------------
function [dataCell, p, stats] = coexpConnUncon (C, G, anatomyPart, dorankedTest)
if nargin < 3
    anatomyPart = 'wholeNervousSystem';
    fprintf(1,'Using whole nervous system BY DEFAULT\n');
end

if nargin < 4
    dorankedTest = true;
    fprintf(1,'Using rank test BY DEFAULT\n');
end

coexpression = GiveMeCoexpression(G,[],true); 
Adj = GiveMeAdj(C,'zeroBinary');


numNeurons = size(Adj,1);

switch anatomyPart
    case 'wholeNervousSystem'
        coexpression = GiveMeCoexpression(G,[],true);  %Pearson_noLR;
    case 'headOnly'
        coexpression = GiveMeCoexpression(G,[],true); 
        head = C.RegionM(:,9);
        Adj = Adj(head == 1, head == 1);
        coexpression = coexpression(head == 1, head == 1);
end

NaNmatrix = nan(numNeurons,numNeurons);
NaNmatrix = tril(NaNmatrix,-1);
Adj = triu(logical(Adj+Adj'),1)+0;
Adj = Adj+NaNmatrix;

% get a mask for unconnected links
Unc = Adj == 0;
% get a mask for connected links
Con = Adj == 1;

Uncon = coexpression.*Unc;
Con = coexpression.*Con;

Uncon(1:size(Uncon,1)+1:end) = nan;
Con(1:size(Con,1)+1:end) = nan;

% exclude NaNs and zeros
Uncon = nonzeros(Uncon);
Uncon=(Uncon(~isnan(Uncon)));
Con = nonzeros(Con);
Con=(Con(~isnan(Con)));

% compare distributions
if dorankedTest==false
    [h,p,ci,stats] = ttest2(Con,Uncon, 'Vartype','unequal');
elseif dorankedTest==true
    [p,h,stats] = ranksum(Con,Uncon);
end

dataCell{1,1} = Con; dataCell{2,1} = Uncon;
% xLabels = {'Connected', 'Unconnected'};
%
% JitteredParallelScatter(dataCell)
%
% set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% set(gca,'box','off');
% ylabel('Coexpression','FontSize', 15);
end
