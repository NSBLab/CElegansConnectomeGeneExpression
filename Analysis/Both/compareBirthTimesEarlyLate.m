% compare coexpression between for early born neurons, late born neurons
% and mixed pairs
connected = true; 
if ~connected
    Adj = ones(279); 
else
Adj = GiveMeAdj(C, 'zeroBinary'); 
Adj = triu(Adj|Adj',1); 
end
maskLR = GiveMeLRMask(C); 
coexpression = GiveMeCoexpression(G,[],false).*~maskLR.*Adj;
BirthTime = C.BirthTime; 

earlyThr = 550;
early = double(BirthTime<earlyThr); 
early = early*early'; 
late = double(BirthTime>earlyThr); 
late = late*late'; 
mix = ~early & ~late; 


coexpEarly = nonzeros(triu(coexpression.*early,1)); 
coexpLate = nonzeros(triu(coexpression.*late,1)); 
coexpMix = nonzeros(triu(coexpression.*mix,1)); 

dataCell{1} = coexpEarly;
dataCell{2} = coexpLate;
dataCell{3} = coexpMix;

[P.earlyLate, ~, S.earlyLate] = ranksum(coexpEarly,coexpLate); 
[P.earlyMix, ~, S.earlyMix] = ranksum(coexpEarly,coexpMix); 
[P.LateMix, ~, S.LateMix] = ranksum(coexpLate,coexpMix); 

JitteredParallelScatter(dataCell); 
xLabels = {'Early', 'Late', 'Mix'};
set(gca,'Xtick', [1 2 3], 'XTickLabel',xLabels, 'FontSize', 15);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 12);
set(gca,'box','off');
ylabel('Gene coexpression, r_\phi','FontSize', 15);

