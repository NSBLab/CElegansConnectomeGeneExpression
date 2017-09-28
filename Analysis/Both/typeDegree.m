function [meanCoexp, numConnections, deg, types, stats] = typeDegree(C,G, types, doMean)
% ------------------------------------------------------------------------------
% Function calculates coexpression for connections involving different
% types of neurons (interneurosn, motor neurons, sensory neurons).
% At each degree threshold we calculate degree coexpression between neurons where
% at least one neuron in a  pair is of a selected type and the degree of both neurons
% is higher than a threshold k.
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% G (gene data structure)
% types (10 13 17) = default (interneurons, motor neurons, sensory neurons)
% coexpMeasure - bydefault is Pearson_noLR
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute 2017-04-10
% ------------------------------------------------------------------------------% ------------------------------------------------------------------------------

if nargin < 3
    types = [10 13 17];% 10 - interneurons, 13 - motor, 17 - sensory;
    fprintf('Interneurons, motor neurons and sensory neurons by DEFAULT\n')
end
if nargin < 4
    doMean = 0;
end

% choose chemical synapse network
Adj = GiveMeAdj(C,'zeroBinary');
coexp = GiveMeCoexpression(G,'',false);
maskLR = GiveMeLRMask(C);
coexp = coexp.*~maskLR; 
fprintf('Using DEFAULT coexpression measure\n');
[~,~,deg] = degrees_dir(Adj);
Adj = triu(Adj|Adj',1); 

% initiate variables
meanCoexp = zeros(max(deg), length(types));
numConnections = zeros(max(deg), length(types));

for t=1:length(types)
    type = C.RegionM(:,types(t));
    % get a mask for type
    B = repmat(type,[1,279]);
    maskType = B|B'; % at least one neuron from a pair needs to be of a selected type
    for d=1:max(deg)
        isHub = (deg>d);
        
        maskDeg = isHub.*isHub'; % both neurons in a pair need to be with degree >k
        % get overall mask
        mask = Adj&maskType&maskDeg; % existing connection & selected type & deg>k
        numConnections(d,t) = sum(mask(:)); % get the number of existing connections for this mask
        coexpMask = coexp.*mask;
        %coexpMask(coexpMask==0)=NaN; % replace zero values with nan
        %coexpMask(isnan(coexpMask)) = [];
        linkDataSpecial = nonzeros(coexpMask(:));
        maskNonspecial = ~mask & Adj;
        linkDataNotSpecial = nonzeros(coexp(maskNonspecial));
        
        if ~isempty(linkDataSpecial(~isnan(linkDataSpecial))) && ~isempty(linkDataNotSpecial(~isnan(linkDataNotSpecial)))
            [p,~,~] = ranksum(linkDataSpecial,linkDataNotSpecial, 'tail', 'right');
            stats(d,t) = p;
        else
            stats(d,t) = NaN;
        end
        
        if doMean
            meanCoexp(d,t) = mean(linkDataSpecial);
        else
            meanCoexp(d,t) = median(linkDataSpecial);
        end
        
    end
end

end
