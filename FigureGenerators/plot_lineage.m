%===============================================================================
% MAIN TEXT:
%===============================================================================
% Plot lineage distance for rich and non-rich links as overlaping
% histograms.
%-------------------------------------------------------------------------------

%%
connected = 1;

D = GiveMeDefault();
kHub = D.kHub;

Adj = GiveMeAdj(C,'zeroBinary');
[~,~,deg] = degrees_dir(Adj);

if ~connected
    Adj=ones(279); 
end

indH = deg>kHub;
maskH = zeros(C.numNeurons,C.numNeurons);

maskH(indH==1, indH==1)=1;
maskP=~maskH;

Adj=Adj|Adj'; 
%Adj = ones(279,279); 
Adj = triu(Adj); 

RelH = C.LineageDistance_noLR(Adj&maskH);
RelH(isnan(RelH))=[];
RelP = C.LineageDistance_noLR(Adj&maskP);
RelP(isnan(RelP))=[];

RelH(isnan(RelH))=0; nH = sum(logical(RelH(:)));
RelP(isnan(RelP))=0; nP = sum(logical(RelP(:)));


[rgb_colorMatrix,labels] = GiveMeColors('RichNONrich');

figure; histogram(RelH, 'Normalization', 'probability', 'LineWidth', 1.5, 'EdgeColor',rgb_colorMatrix(1,:),...
    'FaceColor', brighten(rgb_colorMatrix(1,:),+0.3)); hold on;
histogram(RelP, 'Normalization', 'probability', 'LineWidth', 1.5,'EdgeColor',rgb_colorMatrix(2,:),...
    'FaceColor', brighten(rgb_colorMatrix(2,:),+0.3));
set(gcf,'color','w');box off
set(gca,'Xtick', [0 5 10 15 20 25], 'XTickLabel',[0 5 10 15 20 25], 'FontSize', 15);
set(gca,'Ytick', [0 0.1 0.2 0.3 0.4], 'YTickLabel',[0 0.1 0.2 0.3 0.4], 'FontSize', 15);
xlabel('Lineage distance', 'FontSize', 15); ylabel('Normalised frequency', 'FontSize', 15);
leg = legend(sprintf('Rich (%d pairs)',nH), sprintf('Non-rich (%d pairs)', nP)); legend boxoff ;
set(leg,'FontSize',15);

[ph,~,statsh] = ranksum(RelH,RelP);


