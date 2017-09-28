function [k,sortedLabels] = NeurotransmitterAnal(C,doPlot)
% ------------------------------------------------------------------------------
% Function plots the proportion of neurons with a specific neurotransmitter type as a function of degree.
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% ------------------------------------------------------------------------------
% Ben Fulcher 2017-03-20
% ------------------------------------------------------------------------------

if nargin < 1
    load('CElegansConnectivityData.mat','C');
end
if nargin < 2
    doPlot = true;
end

% Chemical, binary:
Adj = GiveMeAdj(C,'zeroBinary','all');

% Compute degree:
k_in = sum(Adj,1);
k_out = sum(Adj,2)';
k = k_in + k_out;
kmax = max(k);
kRange = min(k):max(k);

% Assign neurotransmitter types:
NTTypeTable = ImportNeurotransmitterTypes();

% Match to regions:
neuronNames = {C.RegionStruct.acronym};
hasZero = find(cellfun(@(x)strcmp(x(3),'0'),neuronNames));
for i = 1:length(hasZero)
    neuronNames{hasZero(i)}(3) = '';
end
% Manually fix DB1/3
fixMe = strcmp(neuronNames,'DB1');
neuronNames{fixMe} = 'DB1/3';
fixMe = strcmp(neuronNames,'DB3');
neuronNames{fixMe} = 'DB3/1';

[~,b,ix_c] = intersect(neuronNames,NTTypeTable.neuronName,'stable');
fprintf(1,'%u/%u neurons match a neurotransmitter label from Pereira (n=%u)\n',...
                length(b),length(neuronNames),height(NTTypeTable));

noMatch = setxor(1:length(neuronNames),b);
% neuronNames(noMatch)

%-------------------------------------------------------------------------------
% Now labels:
sortedLabels = categorical(NTTypeTable.neuroTransmitter(ix_c));
theLabels = categories(sortedLabels);
numLabels = length(theLabels);

%-------------------------------------------------------------------------------
% Plot of the proportions of neurons (set with degree > k) as a function of k
%-------------------------------------------------------------------------------
if doPlot
    f = figure('color','w');
    subplot(4,1,1)
    N = arrayfun(@(x)sum(k==x),kRange);
    histogram(k,100, 'EdgeColor','k','FaceColor','k'); box off;
    xlim([min(k)-0.5,max(k)+0.5]);
    ylabel('Frequency')
    subplot(4,1,2:4)
    % Give colors:
    colors = BF_getcmap('spectral',numLabels,1);
    kRange = min(k):max(k);
    for i = 1:length(kRange)
        isHub = (k > kRange(i));
        categoriesHere = sortedLabels(isHub);
        proportions = countcats(categoriesHere)/length(categoriesHere);
        for j = 1:length(proportions)
            if proportions(j) > 0
                rectangle('Position',[kRange(i),sum(proportions(1:j-1)),1,proportions(j)], ...
                    'FaceColor',colors{j},'EdgeColor',colors{j})
            end
        end
    end
    ylim([0 1]);
    xlabel('Degree, k', 'FontSize', 15);
    ylabel('Proportion of neurons', 'FontSize', 15);
    % Add text labels:
    proportions = countcats(sortedLabels)/length(sortedLabels);
    for j = 1:length(proportions)
        text(kRange(1),sum(proportions(1:j-1))+0.5*proportions(j),theLabels{j})
    end
end

end
