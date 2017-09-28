% do spearman correlation between coexpression and lineage for connected
% neurons
maskLR = GiveMeLRMask(C);
Adj = GiveMeAdj(C,'zeroBinary'); 
Adj = triu(Adj|Adj'); 

coexp = GiveMeCoexpression(G).*~maskLR.*Adj;
Lineage = C.LineageDistance.*Adj.*~maskLR;
figure; scatter(nonzeros(Lineage), nonzeros(coexp));
[r,p]=corr(nonzeros(Lineage), nonzeros(coexp),'type', 'Spearman');