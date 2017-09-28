function D = GiveMeDist(C,whichDistances,excludeLR)

if nargin < 2 || isempty(whichDistances)
    defaults = GiveMeDefault();
    whichDistances = defaults.whatDistance;
end
if nargin < 3
    excludeLR = false;
end
%-------------------------------------------------------------------------------

if excludeLR
    D = C.(sprintf('%s_noLR',whichDistances));
else
    D = C.(whichDistances);
end

end
