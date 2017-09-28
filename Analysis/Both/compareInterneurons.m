% all interneurons
Adj = GiveMeAdj(C,'zeroBinary');
D = GiveMeDefault; 
hub = D.kHub;

coexp = GiveMeCoexpression(G,[],false); 
maskLR = GiveMeLRMask(C); 
coexp = coexp.*~maskLR; 

[~,~,deg] = degrees_dir(Adj); 
Adj = triu(Adj|Adj'); 
inter = C.RegionM(:,10); 
inter = inter|inter'; 
 
isHub = deg>hub; 
isHub = isHub|isHub'; 

special = Adj&isHub&inter; 
nonspecial = Adj&~isHub&inter; 

specialCoexp = coexp(special); 
nonspecialCoexp = coexp(nonspecial); 


[p,~,S] = ranksum(specialCoexp,nonspecialCoexp)
