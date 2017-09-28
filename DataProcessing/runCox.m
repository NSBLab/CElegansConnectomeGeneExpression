
cd ('/projects/kg98/aurina')
load('CelegansGeneDataWS256CertainEmptyEnrichedPartialFINAL.mat'); 

X = G.GeneExpData.Direct; 
coX = CoExpressionPositiveMatch(X,true); 
save('coX', 'coX'); 