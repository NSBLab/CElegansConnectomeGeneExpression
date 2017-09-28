
%===============================================================================
% MAIN TEXT:
%===============================================================================
% First is the topological rich-club curve
%-------------------------------------------------------------------------------
% Parameters:
D = GiveMeDefault();
networkType = 'bd';
nullModel = 'randmio_dir';
numIter = 50;
numNullNetworks = 1000;
doRecompute = false;

% ----(RECOMPUTE)----:
if doRecompute
    networkRC(C,D.whatConns,networkType,nullModel,numIter,numNullNetworks);
end

% ----(LOAD AND PLOT)----
PlotRichClub(D.whatConns,networkType,nullModel,numIter,numNullNetworks);

% Mean distance as a function of a degree
% RichClub(C,G,'connection distance');
% ylim([0 1]);

%% neutortansmitter as a function of degree
NeurotransmitterAnal(C)

%% neuron type as a function of degree
propTypesDegree(C)

%===============================================================================
% SUPP--WEIGHTED RICH CLUBS
%===============================================================================
% For different rich club definitions
networkType = 'wd';
numIter = 50;
numNullNetworks = 1000;

% ----(RECOMPUTE)----:
if doRecompute
    nullModel = 'randmio_dir';
    networkRC(C,D.whatConns,networkType,nullModel,numIter,numNullNetworks);
    nullModel = 'shuffleWeights';
    networkRC(C,D.whatConns,networkType,nullModel,numIter,numNullNetworks);
end

% plot weighted RC (topology fixed)
nullModel = 'shuffleWeights';
h_RCcurveW = PlotRichClubW(D.whatConns,networkType,nullModel, numIter,numNullNetworks, [1 .43 .29; 1 .7 .28], true);
hold on;
% plot mixed RC (topology and weights changed)
nullModel = 'randmio_dir';
h_RCcurveM = PlotRichClubW(D.whatConns,networkType,nullModel, numIter,numNullNetworks, [.6 .73 .45; .91 1 .86],false);

legend([h_RCcurveW,h_RCcurveM],{'\Phi_{normWeighted}','\Phi_{normMixed}'},'Location','NorthWest', 'FontSize', 10);
