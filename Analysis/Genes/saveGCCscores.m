% calculate scores and do FDR correction
% make a list of TOP genes after FDR correction
%% Connected
% function saveGCCscores(chooseType,kHub)
chooseType = 'Connected';
D = GiveMeDefault();
kHub = D.kHub;

[Tcon,~] = GCCbinomial(C,G,chooseType,true, false);
TsortCON = sortrows(Tcon, 'geneScores_corr');
TsortCON.Properties.VariableNames = {'Gene_name' 'pVal' 'Matches_on_link_type' 'All_possible_matches' 'FDR_pVal'};

fileNameOut = fullfile('Data','ermineJdata',sprintf('Correctedscores%s_k%dN10ATLEASTn.txt',chooseType,kHub));
    writetable(TsortCON,fileNameOut,'Delimiter','\t')

%% RichFeeder
chooseType = 'RichFeeder';

[Trf,~] = GCCbinomial(C,G,chooseType,true, false);
TsortRF = sortrows(Trf, 'geneScores_corr');

fileNameOut = fullfile('Data','ermineJdata',sprintf('Correctedscores%s_k%dN10ATLEASTn.txt',chooseType,kHub));
    writetable(TsortRF,fileNameOut,'Delimiter','\t')
