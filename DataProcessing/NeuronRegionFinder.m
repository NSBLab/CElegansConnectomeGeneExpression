function C = NeuronRegionFinder(C,beVocal)
%-------------------------------------------------------------------------------
% 1. Reads in the hierarchy file (from RetrieveHierarchy.py)
% 2. Creates a translation table, and assigns all neurons to specific regions
% 3. Outputs results to matlab file: CElegans_Connectivity_Data.mat
%-------------------------------------------------------------------------------

if nargin < 1
    load('CElegansConnectivityData.mat','C');
end
if nargin < 2
    beVocal = true;
end

% Extract neuron ID
neuronID = cell2mat({C.RegionStruct.id});
numNeurons = length(neuronID);

%-------------------------------------------------------------------------------
%% Create translation table for parents -> children
%-------------------------------------------------------------------------------

% Read in and parse data from hierarchy.csv
[parentName,childName,parentID,childID] = ReadInHierarchy();

% Find anatomy terms annotated directly under the Neuron anatomy term
index = (parentID==3679); % WBbt:0003679
neuronAnatomyIDs = childID(index);
neuronAnatomyNames = childName(index);
numNeuronAnatomyIDs = length(neuronAnatomyIDs);

fprintf(1,'%u direct neuron anatomy terms in hierarchy\n',numNeuronAnatomyIDs);

% Find subsets of anatomy terms to assign labels to each neuron
RegionM = zeros(numNeurons,numNeuronAnatomyIDs);

% Cycle through neurons:
for m = 1:numNeurons
    if beVocal
        fprintf(1,'Neuron %u/%u\n',m,numNeurons);
    end
    % Keep going up the hierarchy, matching, until at the top
    RegionM(m,:) = PropagateUp(neuronID(m));
end
% add information about neuron type according to wormatlas for neurons that
% are not assigned to be neither motor/interneuron/sensory neuron. see: http://www.wormatlas.org/neurons/Individual%20Neurons/PLNframeset.html
% ALA, ALNL, ALNR, AUAL, AUAR - interneurons;
% ALA, CEPDL, CEPDR, CEPVL, CEPVR, PLNL PLNR - sensory neurons
RegionM([23 26 27 52 53],10) = 1;
RegionM([23 83 84 85 86 148 149],17) = 1;

isInter = RegionM(:,strcmp(neuronAnatomyNames,'interneuron'));
isSensory = RegionM(:,strcmp(neuronAnatomyNames,'sensory neuron'));
isMotor = RegionM(:,strcmp(neuronAnatomyNames,'motor neuron'));

multi = isInter+isSensory+isMotor; % check if a neuron belongs to more then one category

RegionM(:,21) = 0; % create empty variable
RegionM(multi>1,21) = 1; % populate with ones if neuron is assigned to more than 1 category.
neuronAnatomyNames{21} = 'multimodal';
% edit information about head-tail assignment according to wormatlas. 
% ALA (23), AVFL(62), AVFR(63), AVG(64), RIFL(173), RIFR(174), RIGL(175),
% RIGR 176), SABD(206), SABVL(207), SABVR(208), SMDVL(225) - in wormatlas
% locations of these neurons are labeled as head. 
% Hierarchy file assignes neither head, nor tail. 
% manualy assign those neurons to head. 
RegionM([23 62 63 64 173 174 175 176 206 207 208 225],9)=1;
C.RegionM = RegionM;
C.neuronAnatomyNames = neuronAnatomyNames;
C.neuronAnatomyIDs = neuronAnatomyIDs;


%-------------------------------------------------------------------------------
% Save data to file: CElegans_Connectivity_Data.mat
% fileName = 'CElegansConnectivityData.mat';
% save(fileName,'C','-append');
% fprintf(1,'Saved neuronAnatomy term annotation results to %s\n',fileName);

%-------------------------------------------------------------------------------
function annotations = PropagateUp(idStart)
    % start at an ID and propagate up the hierarchy, annotating matches on the way
    idHere = idStart;
    annotations = zeros(1,numNeuronAnatomyIDs);

    % find matches to current ID
    if any(idHere==neuronAnatomyIDs)
        annotations(idHere==neuronAnatomyIDs) = 1;
        if beVocal
            fprintf(1,'Annotated for %u\n',idStart);
        end
    else
        % Go up hierarchy to parents
        idUp = parentID(childID==idHere);
        numParents = length(idUp);
        if numParents == 0
            if beVocal
                fprintf(1,'Got to top for %u\n',idStart);
            end
        else
            if beVocal
                fprintf(1,'Going up for %u parents\n',numParents);
            end
            for j = 1:numParents
                % go up hierarchy for each parent:
                annotate_j = PropagateUp(idUp(j));
                annotations(annotate_j==1) = 1;
            end
        end
    end
end

end
