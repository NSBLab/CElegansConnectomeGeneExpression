%===============================================================================
% MAIN TEXT:
%===============================================================================
% Plot modularity-related analysis plots
%-------------------------------------------------------------------------------

miniVersion = true; % Just plot the intra-modular hub distributions
reCompute = false; % recompute module calculation

% Parameters for Louvain algorithm:
gamma = 1;
tau = 0.1;
numRuns = 1000;
onlyConnected = true; % only look at coexpression where connections exist
measureToCompare = 'coexpression';
whatSolution = 'consensus'; % 'best','consensus','ERMM'

% Setting project defaults:
D = GiveMeDefault();
whatMatrix = D.whatConns;
simMatrix = GiveMeAdj(C,'zeroBinary',whatMatrix);
kHub = D.kHub;
%simMatrix = C.Adj_B{3};
%simMatrix = logical(simMatrix+simMatrix')+0;

%-------------------------------------------------------------------------------
% Re-run module calculation and save result to file:
%-------------------------------------------------------------------------------
if reCompute
    [Modules,Q,concModules] = GiveMeModules(C,simMatrix,gamma,tau,numRuns,whatMatrix);
    % Saves to file that is loaded by CoexpressionDistributionsModules
end

%===============================================================================
% Generate a figure
%===============================================================================
% %-------------------------------------------------------------------------------
% % Results for ERMM modules:
% %-------------------------------------------------------------------------------
% CoexpressionDistributionsModules(C,G,'intraInter',onlyConnected,'ERMM',whatMatrix);

f = figure('color','w');
if miniVersion
    ax1 = subplot(1,2,1);
    CoexpressionDistributionsModules(C,G,'hubIntra',onlyConnected,whatSolution,whatMatrix);
    ax1.YTick = 0:0.2:1;
    ax1.YTickLabel = 0:0.2:1;
    box off;
    ax2 = subplot(1,2,2);
    CoexpressionDistributionsModules(C,G,'hubInter',onlyConnected,whatSolution,whatMatrix);
    linkaxes([ax1,ax2])
    ax1.YTick = ax2.YTick;
    ax1.YTickLabel = ax2.YTickLabel;
    f.Position = [168    1025   608   201];
    LabelCurrentAxes('A',ax1,16,'topLeft','k')
    LabelCurrentAxes('B',ax2,16,'topLeft','k')
else
    f.Position = [1000,200,1400,700];

    % Don't plot each point
    extraParams = struct('customSpot','');

    % Plot coexpression distributions within and between modules:
    ax2 = subplot(2,4,3:4);
    CoexpressionDistributionsModules(C,G,'intraInter',onlyConnected,whatSolution,whatMatrix);
    LabelCurrentAxes('B',gca,26,'topRight','k')
    % ax2.Position = [ax2.Position(1),ax2.Position(2)*0.95,ax2.Position(3)*0.5,ax2.Position(4)*0.8]; % [left bottom width height]
    ax2.YTick = 0:0.2:1;
    ax2.YTickLabel = 0:0.2:1;
    ax2.Position = [0.5401    0.53    0.2    0.3];
    box off;

    % Plot coexpression distributions within M3 (hub module) for rich links and
    % non-rich links:
    ax3 = subplot(2,4,7:8);
    CoexpressionDistributionsModules(C,G,'hubIntra',onlyConnected,whatSolution,whatMatrix);
    LabelCurrentAxes('C',gca,26,'topRight','k')
    % ax3.Position = [ax3.Position(1),ax3.Position(2)*2.2,ax3.Position(3)*0.5,ax3.Position(4)*0.8]; % [left bottom width height]
    ax3.YTick = 0:0.2:1;
    ax3.YTickLabel = 0:0.2:1;
    ax3.Position = [0.5401    0.2    0.2    0.3];
    box off;

    % Plot matrix reordered according to modules and coloured for link type:
    ax1 = subplot(2,4,[1,2,5,6]); hold on
    PlotModuleMatrix(C,whatSolution,false);
    LabelCurrentAxes('A',gca,26,'topLeft','w')
    % ax1.Position = [ax1.Position(1)*0.15,ax1.Position(2)*2,ax1.Position(3)*1.8,ax1.Position(4)*1.8]; % [left bottom width height]

    % Do resampling analysis:
    CoexpressionDistributionsModules(C,G,'resample',onlyConnected,whatSolution,whatMatrix);
end
