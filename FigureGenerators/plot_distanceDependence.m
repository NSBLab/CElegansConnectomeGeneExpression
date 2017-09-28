%-------------------------------------------------------------------------------
%% Supp fig on hub-hub links in different anatomical divisions
%-------------------------------------------------------------------------------
% f = figure('color','w'); hold on
% numThresholds = 11;
% doSymmetric = false;
% distinguishClasses = true;
% probabilityVSdistance(C,numThresholds,doSymmetric,distinguishClasses,false);
%-------------------------------------------------------------------------------
hubConnectionDistance(C,'bar');
set(gcf,'Position',[1178,763,456,261]);
% LabelCurrentAxes('C',gca,18,'topLeft');
% set(gca,'XTickLabel',[]); xlabel('')
% subplot(2,1,2)
% justConnected = false;
% makeNewFigure = false;
% whatCorr = 'Pearson_noLR';
% coexpressionDistance(C,G,justConnected,whatCorr,makeNewFigure);
% f.Position = [425    92   557   507];
%-------------------------------------------------------------------------------
