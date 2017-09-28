function D = GiveMeDefault(whatMatrix)

if nargin < 1
    whatMatrix = 'all';
end

D.whatDistance = 'Eucl_Dist_2D';
D.coexpMeasure = 'Pearson';

switch whatMatrix
case 'all'
    D.whatConns = 'all';
    D.kHub = 44;
    fprintf(1,'kHub 44 by default\n');
case 'synaptic'
    D.whatConns = 'ch';
    D.kHub = 41;
    fprintf(1,'kHub 41 by default\n');
end

end
