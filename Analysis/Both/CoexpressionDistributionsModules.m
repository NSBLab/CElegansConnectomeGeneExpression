function CoexpressionDistributionsModules(C,G,plotWhat,onlyConnected,whatSolution,whatMatrix)

if nargin < 3
    plotWhat = 'intraInter';
end
if nargin < 4
    onlyConnected = true;
end
if nargin < 5
    whatSolution = 'consensus';
end
if nargin < 6
    D = GiveMeDefault();
    whatMatrix = D.whatConns;
end

%-------------------------------------------------------------------------------
% Load in data
%-------------------------------------------------------------------------------
gamma = 1;
tau = 0.1;
numRuns = 1000;
D = GiveMeDefault();
kHub = D.kHub;
fileName = sprintf('ModuleAssignment_%smatrix_gamma%.2f_tau%.2f_runs%u.mat',whatMatrix,gamma,tau,numRuns);

switch whatSolution
case 'ERMM'
    moduleAssignment = ERMMFit(C);
    load(fileName,'simMatrix');
    fprintf(1,'Getting module assignment from ERMM\n');
case {'best','consensus'}
    load(fileName)
    fprintf(1,'Loaded data from %s\n',fileName);
    % Get module assignments
    switch whatSolution
    case 'best'
        [~,isBest] = max(Q);
        moduleAssignment = Modules(:,isBest);
        fprintf(1,'Using optimal Q solution\n');
    case 'consensus'
        moduleAssignment = concModules;
        fprintf(1,'Using consensus solution\n');
    end
    % ix = reorder_mod(simMatrix,moduleAssignment);
    % modAssOrd = moduleAssignment(ix);
end

numModules = max(moduleAssignment);
fprintf(1,'%u modules\n',numModules);

%-------------------------------------------------------------------------------
% Characterize hubs across modules
%-------------------------------------------------------------------------------
[~,~,deg] = degrees_dir(simMatrix);
isHub = double(deg > kHub);
for i = 1:numModules
    fprintf(1,'Module %u (%u neurons): %u hubs\n',i,sum(moduleAssignment==i),sum(isHub(moduleAssignment==i)));
end

%-------------------------------------------------------------------------------
% Coexpression:
%-------------------------------------------------------------------------------
coExp = GiveMeCoexpression(G,'',true);
extraParams = struct('customSpot','.');

simMatrixEither = (simMatrix | simMatrix');
simMatrixEither(tril(true(size(simMatrix)))) = 0;

% Make masks
intraModule = repmat(moduleAssignment,[1,C.numNeurons])==repmat(moduleAssignment',[C.numNeurons,1]);
interModule = ~intraModule;
if onlyConnected
    fprintf(1,'Only looking at coexpression for pairs of neurons that are connected\n');
    intraModule(~simMatrixEither) = 0;
    interModule(~simMatrixEither) = 0;
else
    % Look at pairwise:
    fprintf(1,'Pairwise analysis; removing lower diagnoal\n');
    intraModule(tril(true(size(intraModule)))) = 0;
    interModule(tril(true(size(interModule)))) = 0;
end

switch plotWhat
case 'intraInter'
    % Ignore NaNs:
    intraModule(isnan(coExp)) = 0;
    interModule(isnan(coExp)) = 0;
    JitteredParallelScatter({coExp(intraModule),coExp(interModule)},1,1,false,extraParams);
    ax = gca;
    ax.XTick = 1:2;
    ax.XTickLabel = {sprintf('intraModular (%u links)',sum(intraModule(:))),sprintf('interModular (%u links)',sum(interModule(:)))};

    [~,p,~,stats] = ttest2(coExp(intraModule),coExp(interModule),'Vartype','unequal');
    fprintf(1,'%s: t = %.2f, p = %g\n',plotWhat,stats.tstat,p);
    [p,~,stats] = ranksum(coExp(intraModule),coExp(interModule));
    fprintf(1,'%s: z = %.2f, p = %g\n',plotWhat,stats.zval,p);

case 'hubIntra'
    [~,~,deg] = degrees_dir(simMatrix);
    isHub = double(deg > kHub);
    isRich = isHub'*isHub;
    intraModuleRich = intraModule & isRich;
    intraModuleNotRich = intraModule & ~isRich;
    if onlyConnected
        fprintf(1,'Only looking at coexpression where there are connections present\n');
        intraModuleRich(~simMatrixEither) = 0;
        intraModuleNotRich(~simMatrixEither) = 0;
    end
    % Ignore NaNs:
    intraModuleRich(isnan(coExp)) = 0;
    intraModuleNotRich(isnan(coExp)) = 0;
    JitteredParallelScatter({coExp(intraModuleRich),coExp(intraModuleNotRich)},1,1,false,extraParams);
    ax = gca;
    ax.XTick = 1:2;
    ax.XTickLabel = {sprintf('intraModular rich (%u links)',sum(intraModuleRich(:))),sprintf('intraModular non-rich (%u links)',sum(intraModuleNotRich(:)))};

    [~,p,~,stats] = ttest2(coExp(intraModuleRich),coExp(intraModuleNotRich),'Vartype','unequal');
    fprintf(1,'%s: t = %.2f, p = %g\n',plotWhat,stats.tstat,p);
    [p,~,stats] = ranksum(coExp(intraModuleRich),coExp(intraModuleNotRich));
    fprintf(1,'%s: z = %.2f, p = %g\n',plotWhat,stats.zval,p);

case 'hubInter'
    [~,~,deg] = degrees_dir(simMatrix);
    isHub = double(deg > kHub);
    isRich = isHub'*isHub;
    interModuleRich = interModule & isRich;
    interModuleNotRich = interModule & ~isRich;
    if onlyConnected
        fprintf(1,'Only looking at coexpression where there are connections present\n');
        interModuleRich(~simMatrixEither) = 0;
        interModuleNotRich(~simMatrixEither) = 0;
    end
    % Ignore NaNs:
    interModuleRich(isnan(coExp)) = 0;
    interModuleNotRich(isnan(coExp)) = 0;
    JitteredParallelScatter({coExp(interModuleRich),coExp(interModuleNotRich)},1,1,false,extraParams);
    ax = gca;
    ax.XTick = 1:2;
    ax.XTickLabel = {sprintf('interModular rich (%u links)',sum(interModuleRich(:))),sprintf('interModular non-rich (%u links)',sum(interModuleNotRich(:)))};

    [~,p,~,stats] = ttest2(coExp(interModuleRich),coExp(interModuleNotRich),'Vartype','unequal');
    fprintf(1,'%s: t = %.2f, p = %g\n',plotWhat,stats.tstat,p);
    [p,~,stats] = ranksum(coExp(interModuleRich),coExp(interModuleNotRich));
    fprintf(1,'%s: z = %.2f, p = %g\n',plotWhat,stats.zval,p);

case 'resample'
    % Analysis to find expectation on hub-hub coexpression given modular membership
    intraModule = repmat(moduleAssignment,[1,C.numNeurons])==repmat(moduleAssignment',[C.numNeurons,1]);
    interModule = ~intraModule;
    % Ignore NaNs:
    intraModule(isnan(coExp)) = 0;
    interModule(isnan(coExp)) = 0;
    % Remove diagonal:
    intraModule(logical(eye(size(intraModule)))) = 0;

    [~,~,deg] = degrees_dir(simMatrix);
    isHub = double(deg > kHub);
    isRich = isHub'*isHub;
    % Ignore NaNs:
    isRich(isnan(coExp)) = 0;

    if onlyConnected
        fprintf(1,'Only looking at coexpression where there are connections present\n');
        isRich(~simMatrix) = 0;
        intraModule(~simMatrix) = 0;
    end

    richIntra = (isRich & intraModule);
    richInter = (isRich & ~intraModule);
    richAre = [sum(richIntra(:)),sum(richInter(:))];
    numRich = sum(richAre);

    numPerms = 1000;
    meanCoEx = zeros(numPerms,1);
    [intra_i,intra_j] = find(intraModule);
    [inter_i,inter_j] = find(interModule);
    for i = 1:numPerms
        myLinks = zeros(numRich,2);
        rp = randperm(length(intra_i));
        myLinks(1:richAre(1),:) = [intra_i(rp(1:richAre(1))),intra_j(rp(1:richAre(1)))];
        rp = randperm(length(inter_i));
        myLinks(richAre(1)+1:end,:) = [inter_i(rp(1:richAre(2))),inter_j(rp(1:richAre(2)))];
        coExpSample = arrayfun(@(x)coExp(myLinks(x,1),myLinks(x,2)),1:length(myLinks));
        meanCoEx(i) = mean(coExpSample);
    end
    meanRich = mean(coExp(logical(isRich)));
    f = figure('color','w'); hold on;
    histogram(meanCoEx)
    plot(meanRich*ones(2,1),[0,50],'r')
    pValZ = 1 - normcdf(meanRich,mean(meanCoEx),std(meanCoEx));
    fprintf(1,'hubhub: %u intra, %u inter, observed: %g; null: %g +/- %g, p = %g\n',...
                    richAre(1),richAre(2),meanRich,mean(meanCoEx),std(meanCoEx),pValZ);
otherwise
    error('Unknown test %s',plotWhat);
end

ylabel('Gene coexpression, r_\phi');
