function [p, S, dataCell] = coexpELCHRecUnidirUncon(C, G, test, whatCorr)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute
% ------------------------------------------------------------------------------
% script for checking coexpression between neurons, connected with:
% electrical, chemical unidirectional, chemical reciprocal connections and
% unconnected neurons
% Default coexpression measure - Pearson.
% ------------------------------------------------------------------------------

% test - ttest or ranksum

if nargin < 3
    test = 'ranksum';
    fprintf(1,'Comparing distributions with ranksum test BY DEFAULT\n');
end

if nargin < 4
    coexpression = GiveMeCoexpression(G,[],true);
else
    coexpression = GiveMeCoexpression(G, whatCorr,true);
end

NaNmatrix = nan(size(coexpression,1),size(coexpression,1));
NaNmatrix = tril(NaNmatrix);

Electrical = GiveMeAdj(C,'zeroBinary','el');
Chemical = GiveMeAdj(C,'zeroBinary','ch');
Adj  = GiveMeAdj(C);
Adj=Adj|Adj';

% electrical
Electrical = Electrical+NaNmatrix;
% chemical reciprocal
Chemical = Chemical+Chemical';
maskRecip = Chemical==2;
ChemicalRec = maskRecip+NaNmatrix;
% chemical unidirectional
maskUnidir = Chemical==1;
ChemicalUnidir = maskUnidir+NaNmatrix;

Uncon = Adj == 0;
Uncon = Uncon+NaNmatrix;

Electrical = coexpression.*Electrical;
ChemicalRec = coexpression.*ChemicalRec;
ChemicalUnidir = coexpression.*ChemicalUnidir;
Uncon = coexpression.*Uncon;

Electrical = nonzeros(Electrical);
Electrical=(Electrical(~isnan(Electrical)));

ChemicalRec = nonzeros(ChemicalRec);
ChemicalRec=(ChemicalRec(~isnan(ChemicalRec)));

ChemicalUnidir = nonzeros(ChemicalUnidir);
ChemicalUnidir=(ChemicalUnidir(~isnan(ChemicalUnidir)));

Uncon = nonzeros(Uncon);
Uncon=(Uncon(~isnan(Uncon)));


%JitteredParallelScatter(dataCell)
%ylabel('Average coexpression'); legend('Connecteced', 'Unconnected');
switch test
    case 'ttest'
        [~,p.ChRecChUni,~,S.ChRecChUni] = ttest2(ChemicalRec,ChemicalUnidir, 'Vartype','unequal');
        [~,p.ElChUni,~,S.ElChUni] = ttest2(Electrical,ChemicalUnidir, 'Vartype','unequal');
        [~,p.ElChRec,~,S.ElChRec] = ttest2(Electrical,ChemicalRec, 'Vartype','unequal');

    case 'ranksum'
        [p.ChRecChUni,~,S.ChRecChUni] = ranksum(ChemicalRec,ChemicalUnidir);
        [p.ElChUni,~,S.ElChUni] = ranksum(Electrical,ChemicalUnidir);
        [p.ElChRec,~,S.ElChRec] = ranksum(Electrical,ChemicalRec);

        %[p.UnidirUncon,~,S.UnidirUncon] = ranksum(Unidir,Uncon);
        %[p.RecUncon,~,S.RecUncon] = ranksum(Rec,Uncon);

end
dataCell{1,1} = Electrical; dataCell{2,1} = ChemicalRec; dataCell{3,1} = ChemicalUnidir; dataCell{4,1} = Uncon;
end
