function ord = SortStructures(C,criterion)
% Orders nodes by some criterion

if nargin < 1
    criterion = 'type'; % (sensory/motor/inter)
end

isInter = C.RegionM(:,strcmp(C.neuronAnatomyNames,'interneuron'));
isSensory = C.RegionM(:,strcmp(C.neuronAnatomyNames,'sensory neuron'));
isMotor = C.RegionM(:,strcmp(C.neuronAnatomyNames,'motor neuron'));

end
