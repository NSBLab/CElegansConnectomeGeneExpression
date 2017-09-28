function TypeConfoundHubs(doWhat,numPerms)
%-------------------------------------------------------------------------------
% Compare whether high degree neurons are more genetically similar, even amongst
% neurons of a similar set of molecular types
%-------------------------------------------------------------------------------
if nargin < 1
    doWhat = 'neurotransmitter'; % 'anatomy', 'neurotransmitter'
end
if nargin < 2
    numPerms = 1e8; % number of permutations for permutation testing
end
%-------------------------------------------------------------------------------

fprintf(1,'%s with %g permutations\n',doWhat,numPerms);

% Additional parameters:
excludeLR = true;
justConns = false; % whether to just look at connections (otherwise looks at all pairs)
meanOrMedian = 'median';

% Get defaults:
D = GiveMeDefault();

% Get C and G, and paths
LoadAllData;

% Label neurons by neurotransmitter type:
switch doWhat
case 'neurotransmitter'
    [~,ntTypes] = NeurotransmitterAnal(C,false);
    mask = GiveMeMask(C,'NeuroTransIntra');
case 'anatomy'
    ntTypes = LabelNeuronType(C,'anatomy');
    mask = GiveMeMask(C,'AnatomyPairsSymmetric');
end
ntLabels = categories(ntTypes);
numNTLabels = length(ntLabels);

%-------------------------------------------------------------------------------
% Connectivity data:
Adj = GiveMeAdj(C,'zeroBinary');
AdjMaskSym = Adj|Adj';
[~,~,k] = degrees_dir(Adj);
numNeurons = length(Adj);
isHub = (k > D.kHub);

% Get coexpression data:
coExpData = GiveMeCoexpression(G,'',excludeLR);

%-------------------------------------------------------------------------------
% Look at whether distributions are affected by neuron type:
theTypes = fieldnames(mask);
numTypes = length(theTypes);
dataCell = cell(numTypes,1);
justConns = false;
for i = 1:numTypes
    if justConns
        dataCell{i} = coExpData(mask.(theTypes{i}) & triu(true(size(coExpData)),+1) & AdjMaskSym);
    else
        dataCell{i} = coExpData(mask.(theTypes{i}) & triu(true(size(coExpData)),+1));
    end
end
isGood = cellfun(@(x)~all(isnan(x)),dataCell);
dataCell = dataCell(isGood);
theTypes = theTypes(isGood);
[~,ix] = sort(cellfun(@nanmean,dataCell),'descend');
JitteredParallelScatter(dataCell(ix),true,true,true);
ax = gca;
ax.XTick = 1:sum(isGood);
ax.XTickLabel = theTypes(ix);

for i = 1:sum(isGood)
    fprintf(1,'%.2f +/- %.2f [median = %.2f] %s\n',nanmean(dataCell{i}),...
                                nanstd(dataCell{i}),...
                                median(dataCell{i}(~isnan(dataCell{i}))),...
                                theTypes{i});
end

%-------------------------------------------------------------------------------
% Give text info about hubs:
fprintf(1,'%u hubs are:\n',sum(isHub));
hubTypes = ntTypes(isHub);
uTypes = unique(hubTypes);
for i = 1:length(uTypes)
    fprintf(1,'%u(/%u) %s\n',sum(hubTypes==uTypes(i)),sum(ntTypes==uTypes(i)),uTypes(i));
end

%-------------------------------------------------------------------------------
% Get null distribution from random sets of sets of the same neurotransmitter types
%-------------------------------------------------------------------------------
rng('default')
middleCoEx = zeros(numPerms+1,1);
fprintf(1,'%u permutations\n',numPerms);
if justConns
    fprintf(1,'~~only looking at where connections exist~~\n');
end
parfor i = 1:numPerms+1
    if i==1
        theHubs = double(isHub);
    else
        % Take a random sample from the given types
        theHubs = zeros(size(isHub));
        for j = 1:length(uTypes)
            isType = find(ntTypes==uTypes(j));
            theHubs(isType(1:sum(hubTypes==uTypes(j)))) = 1;
            if i==2
                fprintf(1,'%u random neurons of type %s are now hubs!!\n',...
                                sum(hubTypes==uTypes(j)),uTypes(j));
            end
        end
        theHubs = double(isHub(randperm(length(isHub))));
    end
    if justConns
        isRich = theHubs'*theHubs & triu(true(size(Adj))) & AdjMaskSym;
    else
        isRich = theHubs'*theHubs & triu(true(size(Adj)));
    end

    switch meanOrMedian
    case 'mean'
        middleCoEx(i) = nanmean(coExpData(isRich));
    case 'median'
        middleCoEx(i) = nanmedian(coExpData(isRich));
    end
end

pVal = mean(middleCoEx(2:end)>=middleCoEx(1));
fprintf(1,'p = %.3g [on %s]\n',pVal,meanOrMedian);

%-------------------------------------------------------------------------------
save(sprintf('%s-%g-%s.mat',doWhat,numPerms,meanOrMedian));

% f = figure('color','w'); ax = gca; hold on
% histogram(middleCoEx(2:end))
% plot(ones(2,1)*middleCoEx(1),ax.YLim,'r')

%-------------------------------------------------------------------------------
% Extra analyses:
%-------------------------------------------------------------------------------
switch doWhat
case 'anatomy'
    % Also check whether hub-hub increased within head and within tail, and
    % for head-tail pairs
    theHubs = double(isHub);
    if justConns
        isRich = theHubs'*theHubs & triu(true(size(Adj))) & AdjMaskSym;
    else
        isRich = theHubs'*theHubs & triu(true(size(Adj)));
    end
    for i = 1:numTypes
        baseMask = mask.(theTypes{i}) & triu(true(size(coExpData)),+1) & AdjMaskSym;
        a = coExpData(baseMask & isRich);
        b = coExpData(baseMask & ~isRich);
        if ~isempty(a) & ~isempty(b)
            [p,h,stats] = ranksum(a,b);
            fprintf(1,'%.2f +/- %.2f (%u rich); %.2f +/- %.2f (%u non-rich); p = %.3g %s\n',...
                nanmean(a),nanstd(a),length(a),nanmean(b),nanstd(b),length(b),p,theTypes{i});
        end
    end
end
