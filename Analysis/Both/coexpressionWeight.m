% test relationship between connection weight and coexpression
Adj = C.Adj_W{3}; 
coexp = GiveMeCoexpression(G, [], false); 

l = logical(Adj);
mask = GiveMeLRMask(C);

AdjMask = triu(l|l' & ~mask);
A = Adj(AdjMask); 
Coexp = coexp(AdjMask); 

figure; scatter(A, Coexp); xlabel('connection weight'); ylabel('Coexpression'); 

[RHO,PVAL] = corr(A,Coexp); 