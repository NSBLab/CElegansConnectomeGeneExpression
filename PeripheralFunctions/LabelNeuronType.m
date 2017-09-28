function neuronLabel = LabelNeuronType(C,labelWhat)

if nargin < 2
    labelWhat = 'type';
end
%-------------------------------------------------------------------------------

switch labelWhat
case 'type'
    typeAnnot = C.RegionM(:,[find(strcmp(C.neuronAnatomyNames,'interneuron')),...
                            find(strcmp(C.neuronAnatomyNames,'sensory neuron')),...
                            find(strcmp(C.neuronAnatomyNames,'motor neuron'))]);

    % isMulti = C.RegionM(:,strcmp(C.neuronAnatomyNames,'multimodal'));
    numNeurons = 279;
    neuronLabel = cell(numNeurons,1);
    for i = 1:numNeurons
        if sum(typeAnnot(i,:)) > 1
            neuronLabel{i} = 'multi';
        elseif sum(typeAnnot(i,:))==0
            neuronLabel{i} = 'unannotated';
        elseif typeAnnot(i,1)==1
            neuronLabel{i} = 'interneuron';
        elseif typeAnnot(i,2)==1
            neuronLabel{i} = 'sensory';
        elseif typeAnnot(i,3)==1
            neuronLabel{i} = 'motor';
        else
            error('error in elseif??')
        end
    end
    neuronLabel = categorical(neuronLabel);

case 'anatomy'
    isHead = logical(C.RegionM(:,strcmp(C.neuronAnatomyNames,'head neuron')));
    isTail = logical(C.RegionM(:,strcmp(C.neuronAnatomyNames,'tail neuron')));
    isBody = logical(~isHead & ~isTail);
    numNeurons = C.numNeurons;
    neuronLabel = zeros(numNeurons,1);
    neuronLabel(isHead) = 1;
    neuronLabel(isBody) = 2;
    neuronLabel(isTail) = 3;
    neuronLabel = categorical(neuronLabel,1:3,{'head','body','tail'});
end

end
