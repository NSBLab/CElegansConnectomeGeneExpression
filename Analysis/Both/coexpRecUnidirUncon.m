% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute
% ------------------------------------------------------------------------------
% script for checking reciplocally/unidirctionally connected and unconnected link coexpression using
% Pearson coexpression
% ------------------------------------------------------------------------------


function [p, S, dataCell] = coexpRecUnidirUncon (C, G, conType, anatomyPart, test)
% test - ttest or ranksum
if nargin < 3
    Adj = GiveMeAdj(C,'zeroBinary');
end

if nargin < 4
    anatomyPart = 'whole';
end

if nargin < 5
    test = 'ranksum';
    fprintf(1,'Comparing distributions with ranksum test BY DEFAULT\n');
end

if nargin < 6
    coexpression = GiveMeCoexpression(G,[],true); 
end


switch anatomyPart
    case 'whole'
    coexpression = GiveMeCoexpression(G,[],true); 
    case 'head'
        head = C.RegionM(:,9);
        Adj = Adj(head == 1, head == 1);
        coexpression = coexpression(head == 1, head == 1);
end

NaNmatrix = nan(size(coexpression,1),size(coexpression,1));
NaNmatrix = tril(NaNmatrix);
Adj = Adj+Adj';
maskRecip = Adj==2;
maskRecip = maskRecip+NaNmatrix;


maskUnidir = Adj==1;
maskUnidir = maskUnidir+NaNmatrix;



% Adj(Adj == 0) = NaN;
Unc = Adj == 0;
Unc = Unc+NaNmatrix;

Uncon = coexpression.*Unc;
Rec = coexpression.*maskRecip;
Unidir = coexpression.*maskUnidir;


% if only electrical synapses are chosen, use one half of the matrix
% (symetrical)
if strcmp(conType, 'el')
    Uncon = triu(Uncon,1);
    Con = triu(Con,1);
end

Rec = nonzeros(Rec);
Rec=(Rec(~isnan(Rec)));
Unidir = nonzeros(Unidir);
Unidir=(Unidir(~isnan(Unidir)));
Uncon = nonzeros(Uncon);
Uncon=(Uncon(~isnan(Uncon)));



%JitteredParallelScatter(dataCell)
%ylabel('Average coexpression'); legend('Connecteced', 'Unconnected');
switch test
    case 'ttest'
        [~,p.RecUnidir,~,S.RecUnidir] = ttest2(Rec,Unidir, 'Vartype','unequal');
        [~,p.UnidirUncon,~,S.UnidirUncon] = ttest2(Unidir,Uncon, 'Vartype','unequal');
        [~,p.RecUncon,~,S.RecUncon] = ttest2(Rec,Uncon, 'Vartype','unequal');
        dataCell{1,1} = Rec; dataCell{2,1} = Unidir; dataCell{3,1} = Uncon;
    case 'ranksum'
        [p.RecUnidir,~,S.RecUnidir] = ranksum(Rec,Unidir);
        [p.UnidirUncon,~,S.UnidirUncon] = ranksum(Unidir,Uncon);
        [p.RecUncon,~,S.RecUncon] = ranksum(Rec,Uncon);
        dataCell{1,1} = Rec; dataCell{2,1} = Unidir; dataCell{3,1} = Uncon;
end
end
