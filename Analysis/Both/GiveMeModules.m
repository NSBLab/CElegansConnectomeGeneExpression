function [Modules,Q,concModules] = GiveMeModules(C,simMatrix,gamma,tau,numRuns, whatMatrix)
% ------------------------------------------------------------------------------
% function produces modular assignments for a given similarity measure
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure),
% tau - threshold for module detection (0-1); higher values hive mode modules (default 0.1).
% numRuns - how many times to run louvain algorithm                           (default 1000).
% ------------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------------
%
% ------------------------------------------------------------------------------

if nargin < 2
    simMatrix = GiveMeAdj(C,'zeroBinary');
    fprintf(1,'Binary chemical matrix by default\n');
end
if nargin < 3
    gamma = 1; % clasical modularity by default
    fprintf(1,'Using gamma = 1 BY DEFAULT\n');
end
if nargin < 4
    tau = 0.1;
    fprintf(1,'Using tau = 0.1 BY DEFAULT\n');
end
if nargin < 5
    numRuns = 1000;
    fprintf(1,'Using 1000 runs BY DEFAULT\n');
end

numNeurons = size(simMatrix,1);

% ------------------------------------------------------------------------------
% Run Louvain algorithm to get modular assignments
% ------------------------------------------------------------------------------

fprintf(1,'Across %u runs of Louvain community detection...',numRuns);
% Assign modular assignment and modularity scores:
Modules = zeros(numNeurons,numRuns);
Q = zeros(numRuns,1);
for jj=1:numRuns
    [M, Q1] = community_louvain(simMatrix, gamma);
    Modules(:,jj) = M;
    Q(jj) = Q1;
end
fprintf(1,' Done\n');

% Modular assignment agreement, weighted by modularity scores:
fprintf(1,'Getting modular assignment agreement...');
D = agreement_weighted(Modules,Q);
fprintf(1,' Done.\n');

%% Consensus assignment of neurons to modules:
fprintf(1,'Getting consensus assignment...');
concModules = consensus_und(D,tau,numRuns);
fprintf(1,' Done.\n');

%===============================================================================
% Save results:
%===============================================================================
fileName = sprintf('ModuleAssignment_%smatrix_gamma%.2f_tau%.2f_runs%u.mat',whatMatrix,gamma,tau,numRuns);
fileName = fullfile('Data',fileName);
save(fileName)
fprintf(1,'Saved module analysis to %s\n',fileName);

end
