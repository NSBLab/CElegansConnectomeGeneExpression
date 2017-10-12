function [p, stats, dataCell] = plotInterneuronsDegree(C,G, coexpMeasure, doSubset, doConnected, kHub)
% ------------------------------------------------------------------------------
% Function plots coexpression between pairs of hub interneurons and non-hub interneurons as a function of degree;
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% G (gene data structure)
% doSubset -
%   true (choose a subset of non-hub interneurons that have similar anatomy to hub interneurons)
%   false (choose all interneurons)
% doConnected -
%   true (compare coexpression for only connected neurons
%   false (compare coexpression for all possible pairs of neurons)
% kHub - hub threshold          default - 41.
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute 2017-03-20
% ------------------------------------------------------------------------------
%%

if nargin < 3 || isempty(coexpMeasure)
   coexpression = GiveMeCoexpression(G,'',true);
         fprintf(1,'Using Pearson coexpression BY DEFAULT\n');
else
   coexpression = GiveMeCoexpression(G,coexpMeasure,true);
end

if nargin < 4
    doSubset = false;
     fprintf('Using all non-hub interneurons by DEFAULT\n')
end

if nargin < 5
    doConnected = true;
     fprintf('Using all links by DEFAULT\n')
end

if nargin < 6
    D = GiveMeDefault();
    kHub = D.kHub;
end
%-------------------------------------------------------------------------------

Adj = GiveMeAdj(C,'zeroBinary');
[~,~,deg] = degrees_dir(Adj);

hubind = find(deg>kHub);
if doSubset == true
    % get nonhub interneurons that have same anatomical properties
    NonHubinterneurons = {'AVFL', 'AVFR', 'AVHL', 'AVHR','AVJL', 'AVJR','AVKR'}; 
    [~,Nhubind] = intersect(C.RegionAcronyms, NonHubinterneurons);
else
    interneurons = find(C.RegionM(:,10));
    Nhubind = setdiff(interneurons,hubind);
    NonHubinterneurons = C.RegionAcronyms(Nhubind);
end

if doConnected
    % take all pairs of connected neurons
    Adj = Adj|Adj';
    coexpression = coexpression.*Adj; %Pearson_noLR.*Adj;
end


hubcoexp = nonzeros(triu(coexpression(hubind, hubind),1));
dataCell{1} = (hubcoexp(~isnan(hubcoexp)));

nonhubcoexp = nonzeros(triu(coexpression(Nhubind, Nhubind),1));
dataCell{2} = (nonhubcoexp(~isnan(nonhubcoexp)));

[p,~,stats] = ranksum(dataCell{1},dataCell{2});
%[~,p,~,stats] = ttest2(dataCell{1},dataCell{2}, 'Vartype', 'unequal');

% xLabels = {'Hub interneurons', 'Non-hub interneurons'};
% JitteredParallelScatter(dataCell)
% set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% set(gca,'box','off');
% ylabel('Coexpression','FontSize', 15);

degH = deg(hubind);
degNH = deg(Nhubind);
Deg = horzcat(degH, degNH);
inds = horzcat(hubind, Nhubind');
mix = (vertcat(Deg, inds))';
smix = sortrows(mix,-1);

if doConnected
    rawcoexp = coexpression.*Adj;
    coexp = rawcoexp(smix(:,2),smix(:,2));
else
    rawcoexp = coexpression;
    coexp = rawcoexp(smix(:,2),smix(:,2));
end
%  figure; imagesc(coexp);
%  colormap([flipud(BF_getcmap('blues',9));[1 1 1]; BF_getcmap('reds',9)]);
%  Names = C.RegionAcronyms(smix(:,2));
%  set(gca, 'XTick', 1:length(Names), 'XTickLabels',Names, 'YTick', 1:length(Names), 'YTickLabels', Names);
%  caxis([-1,1]);colorbar; xtickangle(45); set(gcf,'color','w');
%  title('Coexpression interneurons');

end
