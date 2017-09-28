% Components of main text figure on spatial effects (connection probability
% and gene expression)
f = figure('color','w'); hold on
f.Position = [1000        1149         614         189];
numThresholds = 11;
doSymmetric = false;
distinguishClasses = true;
probabilityVSdistance(C,numThresholds,doSymmetric,distinguishClasses,false);
% SaveAllFigures('SpatialEffects','eps',false);

%-------------------------------------------------------------------------------
% coexpressionDistance(C,G);
