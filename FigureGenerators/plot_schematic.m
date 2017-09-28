%===============================================================================
% SCHEMATIC
%===============================================================================

%f = figure('color','w');

% Connectivity plot:
%subplot(1,4,1)

sortNeurons = 'xPos';
sortGenes = 'descend';
addMetadata = true;
PlotMatrixRichFeederPeripheral(C,[],'',true,sortNeurons,addMetadata);

% Gene expression plot:
 %subplot(1,4,2:3)
f = figure('color','w');
plotGeneExpression(G,C,false,sortGenes,false);

% Coexpression plot
%subplot(1,4,4)
% plotCoexpression(G,C,'Pearson',true,'ugly',sortNeurons);
