function coExp = GiveMeCoexpression(G,whatCorr,excludeLR)

% SET DEFAULTS:
if nargin < 2 || isempty(whatCorr)
    D = GiveMeDefault();
    whatCorr = D.coexpMeasure;
    fprintf(1,'Using %s by default\n',whatCorr);
end
if nargin < 3
    excludeLR = false;
end

%-------------------------------------------------------------------------------
if excludeLR
    coExp = G.Corr.(sprintf('%s_noLR',whatCorr));
    fprintf(1,'Excluding bilateral pairs\n');
else
    coExp = G.Corr.(whatCorr);
end

%-------------------------------------------------------------------------------
% Remove diagonal:
warning('Removing diagonal')
coExp(logical(eye(size(coExp,1)))) = 0;

end
