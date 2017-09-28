
%-------------------------------------------------------------------------------
% Compute binary degree:
Adj = GiveMeAdj(C,'zeroBinary');
[~,~,deg] = degrees_dir(Adj);

%-------------------------------------------------------------------------------
% Compute annotation proportion
gData = [G.GeneExpData.Direct];
annProp = mean(gData,2);

%-------------------------------------------------------------------------------
% Plot:
f = figure('color','w');
plot(deg,annProp,'.k')
xlabel('Degree, k')
ylabel('Proportion of genes annotated')
