function [dataCell,hubModCell,pwithinbetween,zwithinbetween,pwithinbetweenNOR,statswithinbetweenNOR,BTWstats,AdjORD] = ...
                analyzeModulesCoexpression(C,G,gamma,coExpMeasure,tau,numRuns,kHub,measureToCompare,whatSolution)
% ------------------------------------------------------------------------------
% function produces modular assignments and compares coexpression between
% modules as well as within modules vs between modules.
% it also plots:
% number of neurons in each module
% number of hubs in each module
% connectivity matric colored by rich/feeder/peripheral links (original and
% ordered by modules)
% coexpression idstributions within each module
% the proportion of neuron types in each module
% density of neurons in space (head-tail) for each module - spatial distribution.
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure),
% G (gene data structure);
% coExpMeasure - choose coexpression measure as defined in G.Corr. default Pearson_noLR
% tau - threshold for module detection (0-1); higher values hive mode modules (default 0.1).
% numRuns - how many times to run louvain algorithm                           (default 1000).
% kHub - threshold for degining hubs (degree where hubs defined)              (default 42).
% measureToCompare - choose the measure you want to compare between/within modules:
% 'coexpression' or 'lineage';
% ------------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------------
% overlap - how many neurons of each category fall in each module(rows-neuron type, column - modules)
% BTWstats - coexpression difference statistics between pairs of modules
% WOstats - coexpression difference between "within module" to the outside
% of the module
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2017-03-20
% ------------------------------------------------------------------------------

if nargin < 3
    gamma = 1; % clasical modularity by default
    fprintf(1,'Using gamma = 1 BY DEFAULT\n');
end
if nargin < 4
    D = GiveMeDefault();
    coExpMeasure = D.coexpMeasure;
    fprintf(1,'Using Pearson coexpression (excluding L-R matches) BY DEFAULT\n');
end
if nargin < 5
    tau = 0.1;
    fprintf(1,'Using tau = 0.1 BY DEFAULT\n');
end
if nargin < 6
    numRuns = 1000;
    fprintf(1,'Using 1000 runs BY DEFAULT\n');
end
if nargin < 7
    kHub = GiveMeDefault('kHub');
    fprintf(1,'Using hub threshold, k > 41 BY DEFAULT\n');
end
if nargin < 8
    measureToCompare = 'coexpression';
end
if nargin < 9
    whatSolution = 'best';
    whatSolution = 'concensus';
end


%-------------------------------------------------------------------------------
% Load in data saved by GiveMeModules:
%-------------------------------------------------------------------------------
fileName = sprintf('ModuleAssignment_gamma%.2f_tau%.2f_runs%u.mat',gamma,tau,numRuns);
load(fileName)
fprintf(1,'Loaded data from %s\n',fileName);
Adj = simMatrix;

[~,~,deg] = degrees_dir(Adj);
isHub = find(deg>kHub);

switch whatSolution
case 'best'
    [~,isBest] = max(Q);
    moduleAssignment = Modules(:,isBest);
case 'concensus'
    moduleAssignment = concModules;
end

%-------------------------------------------------------------------------------
% reorder connectivity matrix according to modules)
%-------------------------------------------------------------------------------
[ix_mod,~] = reorder_mod(Adj,moduleAssignment);

%-------------------------------------------------------------------------------
%% do analysis
%-------------------------------------------------------------------------------
% plot how manu nodes are in each module
numMod = length(unique(moduleAssignment));
modSize = zeros(numMod,1);
for y = 1:numMod
    modSize(y) = length(find(moduleAssignment==y));
end
% figure; bar(modSize, 'FaceColor',[1 .46 .22],'EdgeColor',[.32 .28 .31],'LineWidth',3); ylabel('Number of neurons', 'FontSize', 15); xlabel('Module', 'FontSize', 15);
% box off;

% get coexpression (or lineage) data
switch measureToCompare
    case 'coexpression'
        measuretest = GiveMeCoexpression(G,coExpMeasure,true); %Pearson_noLR;
        plotName = 'Coexpression';
    case 'lineage'
        measuretest = C.LineageDistance;
        plotName = 'Lineage';
    case 'distance'
        measuretest = GiveMeDist(C);
        plotName = 'Connectivity cost';
end

%-------------------------------------------------------------------------------
% make symetrical connectivity matrix and take the upper triangle to get
% all connections in the connectome (don't care if they resiprocal or not).
Adj = double(Adj+Adj');
Adj(Adj==0)=NaN;  % replace nonexisting connections with NaN
nans = tril(nan(numNeurons,numNeurons)); % make nan matrix and use that to exclude the lower triangle.
Adj = Adj+nans;

measuretest = measuretest.*Adj;

% ------------------------------------------------------------------------------
% TEST difference in coexpression within modules outside modules
% ------------------------------------------------------------------------------

mIn = zeros(numNeurons,numNeurons);
for m1=1:numMod
    for o=1:numNeurons
        for oo=o+1:numNeurons
            if moduleAssignment(o)==m1 && moduleAssignment(oo)==m1
                mIn(o,oo) = measuretest(o,oo);
            end
        end
    end
end

%replace NaNs with 0
mInmask = mIn;
mInmask(isnan(mInmask))=0;
mask = logical(mInmask);
mask = Adj-mask;
mOut = measuretest.*mask;
within = nonzeros(mIn(~isnan(mIn)));
between = nonzeros(mOut(~isnan(mOut)));

% compare distributions
[pwithinbetween,~,stats] = ranksum(within,between);
dataCell{1} = within;
dataCell{2} = between;
% JitteredParallelScatter(dataCell);
% set(gca,'Xtick', 1:2, 'XTickLabel',{'Within modules', 'Between modules'}, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% ylabel('Coexpression', 'FontSize', 15); box off;
zwithinbetween = stats.zval;

% ------------------------------------------------------------------------------
% 2. TEST difference between modules
% ------------------------------------------------------------------------------
%% get coexpression values within each module and coexpression values between each module and other modules
% save them as lists
measureList = cell(numMod,2);
for mod=1:numMod
    % geto coexpresion values within modules
    module = moduleAssignment==mod;
    measureMOD = measuretest(module==1,module==1);
    measureMOD = measureMOD(~isnan(measureMOD));
    measureList{mod,1} = measureMOD;
end

% compare coexpression between modules without rich links in them
hubModCoexpeach = zeros(279,279,numMod);
hubModCoexpeachRICH = zeros(279,279,numMod);

for hubMod = 1:numMod
    for p=1:numNeurons
        for k=1:numNeurons
            if moduleAssignment(p)==hubMod && moduleAssignment(k)==hubMod
                hubModCoexpeach(p,k,hubMod) = measuretest(p, k);
            end
        end
    end
end

for hubMod2 = 1:numMod
    for q=1:numNeurons
        for w=1:numNeurons
            if ~isempty(intersect(isHub,q)) && ~isempty(intersect(isHub,w))
                    hubModCoexpeachRICH(q,w,hubMod2) = hubModCoexpeach(q,w, hubMod2);
            else
                hubModCoexpeachRICH(q,w,hubMod2) = 0;
            end
        end
    end
    modNoR = squeeze(hubModCoexpeach(:,:,hubMod2)) - squeeze(hubModCoexpeachRICH(:,:,hubMod2)) ;
    measureList{hubMod2,2} = nonzeros(modNoR(~isnan(modNoR)));
end
% JitteredParallelScatter(measureList(:,2)); title ('Without Rich links');
% set(gca,'Xtick', 1:numMod, 'XTickLabel',1:numMod, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% ylabel('Coexpression', 'FontSize', 15); xlabel('Modules', 'FontSize', 15);  box off;


% compare coexpression values between pairs of modules using test for
% non-normal distributions (coexpression distributions are non-normal)
pBTWmodall = zeros(numMod,numMod);
zBTWmodall = zeros(numMod,numMod);
pBTWmodNOR = zeros(numMod,numMod);
zBTWmodNOR = zeros(numMod,numMod);
for ll=1:numMod
    for kk=1:numMod
        [p1,~,stats1] = ranksum(measureList{ll,1}, measureList{kk,1}); % test for non-normal distributions.
        [p2,~,stats2] = ranksum(measureList{ll,2}, measureList{kk,2}); % test for non-normal distributions.
        if length(measureList{ll,1})<=20 || length(measureList{ll,1})<=20 % there are no z scores if distributions are smaller than 20 items
            zBTWmodall(ll,kk) = NaN;
            zBTWmodNOR(ll,kk) = NaN;
        else
           % if p1<0.05
                pBTWmodall(ll,kk) = p1;
                zBTWmodall(ll,kk) = stats1.zval;
                pBTWmodNOR(ll,kk) = p2;
                zBTWmodNOR(ll,kk) = stats2.zval;
            %else
%                 pBTWmod(ll,kk) = NaN;
%                 zBTWmod(ll,kk) = NaN;
            %end
        end
    end
end
BTWstats = table(pBTWmodall,zBTWmodall,pBTWmodNOR,zBTWmodNOR);

% ------------------------------------------------------------------------------
% 3. TEST overrepresentation of neuron types in modules
% ------------------------------------------------------------------------------
%% is any group (interneurons/motor/sensory) overrepresented in one of the modules
% check how many neurons fall into each of modules according to a hierarchical group
types =[10 13 17 21];
% 10 - interneurons
% 13 - motor neurons
% 17 - sensory neurons
overlap = zeros(length(types), numMod);

for k=1:numMod
    one = find(moduleAssignment==k);
    for i=1:length(types)
        type = types(i);
        I = find(C.RegionM(:,type));
        overlap(i,k) = length(intersect(I,one));  % number of neurons from each category in each module
    end

end

NrinMOD = sum(overlap,1);
for oo=1:k
   overlap(:,oo) = overlap(:,oo)/NrinMOD(oo);
end

PP1 = {'Interneurons','Motor','Sensory','Multimodal'};
%     figure; h1=bar(overlap', 'stacked'); title('Module composition according to neuron types') ;
%     xlabel('Modules', 'FontSize', 15); ylabel('Proportion', 'FontSize', 15);
%     c = [1 .46 .22; .32 .28 .31; 0 .5 .5; .98 .85 .37];
%     h1(1).Parent.Parent.Colormap = c;
%     legend(h1,PP1);
% ------------------------------------------------------------------------------
% 4. PLOT (no calculations made here)
% ------------------------------------------------------------------------------
% label nodes in connectivity according to module

A = PlotMatrixRichFeederPeripheral(C,kHub,'ch',false,false);
% cout proportion of rich/feeder links in each module
% proportion - what proportion of R\F\P links fall into this module
RFPproportions = zeros(3,numMod);
for tt=1:numMod
    MaskMod = zeros(numNeurons,numNeurons);
    MaskMod(moduleAssignment==tt, moduleAssignment==tt)=1;
    RFPmask = MaskMod.*A;
    MaskMod = logical(RFPmask);
    propR = nansum(nansum(RFPmask==3))/nansum(MaskMod(:)); %nansum(nansum(A==3));
    propF = nansum(nansum(RFPmask==2))/nansum(MaskMod(:)); %nansum(nansum(A==2));
    propP = nansum(nansum(RFPmask==1))/nansum(MaskMod(:)); %nansum(nansum(A==1))
    RFPproportions(:,tt) = [propR, propF, propP];
    propR1 = nansum(nansum(RFPmask==3))/nansum(nansum(A==3));
    propF1 = nansum(nansum(RFPmask==2))/nansum(nansum(A==2));
    propP1 = nansum(nansum(RFPmask==1))/nansum(nansum(A==1));
    RFPproportions1(:,tt) = [propR1, propF1, propP1];
    % how manu hubs in the module
    numHubs(tt) = length(intersect(find(moduleAssignment==tt),isHub));
end
   % plot the number fo hubs in each module
%    figure; bar(numHubs, 'FaceColor', [1 .46 .22],'EdgeColor',[.32 .28 .31],'LineWidth',3); xlabel('Modules', 'FontSize', 15); ylabel('Number of hubs', 'FontSize', 15);
%    set(gca,'Xtick', 1:numMod, 'XTickLabel',1:numMod); box off;

% comare coexpression within module with max number of hubs
% 1. all links
[~, hubMod] = max(numHubs);
hubModCoexpALL = zeros(279,279);


for p=1:numNeurons
    for k=1:numNeurons
        if moduleAssignment(p)==hubMod && moduleAssignment(k)==hubMod
            hubModCoexpALL(p,k) = measuretest(p, k);
        end
    end
end

% 2. only rich links within this module
hubModCoexpRICH = zeros(numNeurons,numNeurons);
for p=1:numNeurons
    for k=1:numNeurons
        if ~isempty(intersect(isHub,p)) && ~isempty(intersect(isHub,k))
            hubModCoexpRICH(p,k) = hubModCoexpALL(p, k);
        end
    end
end
% 3. links without rich links in this module
hubModCoexpNONRICH = hubModCoexpALL - hubModCoexpRICH;
%hubModCell{2} = nonzeros(hubModCoexpALL(~isnan(hubModCoexpALL)));
hubModCell{1} = nonzeros(hubModCoexpRICH(~isnan(hubModCoexpRICH)));
hubModCell{2} = nonzeros(hubModCoexpNONRICH(~isnan(hubModCoexpNONRICH)));
% compare distributions
%[pmodRNR,~,statsmodRNR] = ranksum(hubModCell{1}, hubModCell{3});% rich vs non-rich
[pwithinbetweenNOR,~,statswithinbetweenNOR] = ranksum(hubModCell{1}, hubModCell{2});% all vs non-rich

% JitteredParallelScatter(hubModCell);
% set(gca,'Xtick', 1:3, 'XTickLabel',{'Rich links within module', 'All links within module', 'Non-rich links within module'}, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% ylabel('Coexpression', 'FontSize', 15); box off;

AdjORD = A(ix_mod,ix_mod);
% % plot matrix ordered according to modules and coloured according to
% % rich/feeder/peripheral
% figure('color','w'); imagesc(AdjORD); axis('square'); title('Connectome ordered in modules original directed')
% colors = GiveMeColors('richFeedInOutPeripheral');
% colormap(colors);
% set(gca,'Xtick', [], 'XTickLabel',[]);
% set(gca,'Ytick', [], 'YTickLabel',[]);

% within module
% datawithin = measureList(:,1);
% JitteredParallelScatter(datawithin);
% set(gca,'Xtick', 1:numMod, 'XTickLabel',1:numMod, 'FontSize', 15);
% set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 15);
% ylabel('Coexpression', 'FontSize', 15); xlabel('Modules', 'FontSize', 15);  box off;
%% compare neuron density in space for each module

coordinates = C.Pos(:,1);
NumBins = 100;
edges = linspace(min(coordinates),max(coordinates), NumBins);
densSpace = zeros(numMod,10);
for qq=1:numMod
    module = moduleAssignment==qq;
    neurons = coordinates(module==1);
    Y = discretize(neurons,edges);
    for bin = min(Y):max(Y)
        densSpace(qq,bin) = length(find(Y==bin))/length(Y);
    end
end

% figure; imagesc(fliplr(densSpace)); colormap([flipud(BF_getcmap('blues',9)); BF_getcmap('reds',9)]); colorbar;
%         set(gca, 'XTick', [NumBins*0.1 NumBins*0.9], 'XTickLabels', {'Head', 'Tail'}, 'FontSize', 15);
%         ylabel('Modules', 'FontSize', 15);
%         set(gca,'Ytick', 1:numMod, 'YTickLabel',1:numMod);

end
