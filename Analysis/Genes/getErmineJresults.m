% compile ermineJ results and keep only signifficant values after
% correction
type = 'ORA_RF_FDR0.0001'; 
fileINname = sprintf('%s.erminej.txt', type);  
signThr = 0.05; 

ermineJResultsBPcon = ReadInErmineJ(fileINname); 
isSign = ermineJResultsBPcon.corr_pval<signThr; 
%ermineJResultsBPcon = ermineJResultsBPcon(1:15,:);
% 
% ermineJResultsMFcon = ReadInErmineJ('CON_MF.erminej.txt'); 
% isSign = ermineJResultsMFcon.pVal_corr<signThr; 
% ermineJResultsMFcon = ermineJResultsMFcon(isSign,:);
% 
% ermineJResultsCCcon = ReadInErmineJ('CON_CC.erminej.txt'); 
% isSign = ermineJResultsCCcon.pVal_corr<signThr; 
% ermineJResultsCCcon = ermineJResultsCCcon(isSign,:);

CONall = vertcat(ermineJResultsBPcon); %, ermineJResultsMFcon, ermineJResultsCCcon); 
CONall.corr_pval = mafdr(CONall.pval,'BHFDR',true);
CONall = CONall(1:15,:);
GOtableCON = table(); 
GOtableCON.GOcategory = CONall.GOID;
GOtableCON.Description = CONall.GOName;
GOtableCON.NumGenes = CONall.numGenes; 
num_dig = 4;
GOtableCON.Pval = round(CONall.pval*(10^num_dig))/(10^num_dig);
GOtableCON.Pval_corr = round(CONall.corr_pval*(10^num_dig))/(10^num_dig);

cd('/Users/Aurina/GoogleDrive/Genetics_connectome/CElegans/CElegansCode/Data/ermineJdata/ORAanalysis/ForTables'); 
fileOUTname = sprintf('ermineJresults%s.csv', type); 
writetable(GOtableCON,fileOUTname)


