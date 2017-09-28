function G = Gene_Correlation(G,whatData,annotationType)
%annotationType
%% Taqi Ali 23-7-15
% Takes the genetic expression of each neuron and finds the pearson correlation
% The pearson correlation formula is from Kaufman et al. 'Protocol_S1', and
% applies to binary data
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2016-10-26
% Opition for choosing correlation method added
% (choose: 'Pearson','SimpleMatching', 'YuleQ', 'YuleY', 'Jaccard', 'Dice')
% ------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% Extract Direct annotation matrices
EM_Direct = G.GeneExpData.Direct;
% remove columns with all ones
EM_Direct = EM_Direct(:,~all(EM_Direct,1));
numNeurons = size(EM_Direct,1);

% calculate (expression similarity matrix)
corrMeas = {'Loevinger','YuleQ','YuleY','Pearson','SimpleMatching','Jaccard','Prop','Chi'};
corrMeasLR = {'Pearson','coX'};
numCorrMeas = length(corrMeas);
numCorrMeasLR = length(corrMeasLR);
% Initialize:
for i = 1:numCorrMeas
    G.Corr.(corrMeas{i}) = zeros(numNeurons,numNeurons);
end

% Only loop over upper diagonal -- all definitions should be symmetric
fprintf('Computing pairwise correlations between all pairs of %u neurons...',numNeurons)
for i = 1:numNeurons
    for j = i:numNeurons
        for k = 1:numCorrMeas
            G.Corr.(corrMeas{k})(i,j) = BinaryCorr(EM_Direct(i,:),EM_Direct(j,:),corrMeas{k});
        end
    end
end

% Symmetrize:
for k = 1:numCorrMeas
    lowerTriangle = G.Corr.(corrMeas{k})';
    lowerTriangle(triu(true(numNeurons))) = 0; % ensure diagonal, and upper triangle is zero
    G.Corr.(corrMeas{k})(tril(true(numNeurons),-1)) = 0; % ensure lower diagonal is zero
    G.Corr.(corrMeas{k}) = G.Corr.(corrMeas{k}) + lowerTriangle;
end
fprintf(1,' Done.\n');

%===============================================================================
% ADD additional correlation measure based on possitive matches:
fprintf(1,'Computing positive, match-based correlations using binomial probabilities...');
G.Corr.coX = CoExpressionPositiveMatch(EM_Direct);
fprintf(1,' Done\n');

%===============================================================================
% Create additional matrices with no LR coexpression for each correlation measure
%===============================================================================
maskIsLR = GiveMeLRMask(G);

% Put NaNs over LR symmetric connections:
for p = 1:numCorrMeasLR
    fieldName = sprintf('%s_noLR',corrMeasLR{p});
    G.Corr.(fieldName) = G.Corr.(corrMeasLR{p});
    G.Corr.(fieldName)(maskIsLR) = NaN;

    fprintf(1,'Converted %u entries to NaNs in %s\n',sum(maskIsLR(:)),fieldName);
end

fileNameSave = sprintf('CelegansGeneDataWS%d%s.mat',whatData,annotationType);
fileNameSave = fullfile('Data',fileNameSave);
save(fileNameSave,'G');

fprintf(1,'Direct annotation gene coexpression data saved to %s\n',fileNameSave);

%-------------------------------------------------------------------------------
% Append to the gene data file
% fileName = sprintf('CelegansGeneDataWS%d%s.mat',whatData,annotationType);
% fileName = fullfile('Data',fileName); % place in the Data directory
% save(fileName,'G');
% fprintf(1,'Direct annotation gene coexpression data saved to %s\n',fileName);

%-------------------------------------------------------------------------------
% ---SUBFUNCTIONS---
%-------------------------------------------------------------------------------
    function rho = BinaryCorr(x,y,CorrMethod)
        % 1. Calculate contingency table
        a = sum(x==1 & y==1); % positive matches
        b = sum(x==1 & y==0); % mismatch
        c = sum(x==0 & y==1); % mismatch
        d = sum(x==0 & y==0); % negative matches
        n = a + b + c + d; % total number of observations
        p1 = a+b;
        q1 = c+d;
        p2 = a+c;
        q2 = b+d;
        if n~=length(x)
            error('error computing contingency table');
        end

        % 2. Calculate correlation coefficient
        switch CorrMethod
            case 'Loevinger'
                rho = (a*d-c*b)/(min(p1*q2, p2*q1));
            case 'Pearson'
                rho = (a*d-c*b)/sqrt((a+b)*(c+d)*(a+c)*(b+d));
            case 'SimpleMatching'
                rho = (a+d)/(a+b+c+d);
            case 'YuleQ'
                rho = ((a*d)-(b*c))/((a*d)+(b*c));
            case 'YuleY'
                rho = (sqrt(a*d)-sqrt(b*c))/(sqrt(a*d)+sqrt(b*c));
            case 'Jaccard'
                 rho = a/(a+b+c);
            case 'Dice'
                rho = (2*a)/((2*a)+b+c);
            case 'Prop'
                rho = a/min((a+b),(a+c));
            case 'Chi'
               rho = sqrt((n*(a*d - c*b)^2)/((a+b)*(c+d)*(a+c)*(b+d)));
        end
    end

end
