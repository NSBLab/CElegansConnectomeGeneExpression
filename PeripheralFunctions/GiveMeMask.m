function mask = GiveMeMask(C,whatLinks,networkType,linkedAdj,extraParam,shutUp)
%-------------------------------------------------------------------------------
% Returns a mask of different (labeled) types of links in a structure
%-------------------------------------------------------------------------------

if nargin < 6
    shutUp = true;
end

% ------------------------------------------------------------------------------
if exist('linkedAdj','var')
    if ~islogical(linkedAdj)
        linkedAdj = logical(linkedAdj);
    end
end

switch whatLinks
case 'NeuroTransIntra'
    % Within different neurotransmitterTypes
    [~,ntTypes] = NeurotransmitterAnal(C,false);
    ntLabels = categories(ntTypes);
    numNTLabels = length(ntLabels);
    for i = 1:numNTLabels
        isType = double((ntTypes==ntLabels(i)));
        switch ntLabels{i}
        case 'ACh & 5HTACh & 5HT'
            ntLabels{i} = 'ACh5HTACh5HT';
        case 'ACh & 5HT'
            ntLabels{i} = 'ACh5HT';
        case 'Glu & 5HTGlu & 5HT'
            ntLabels{i} = 'Glu5HTGlu5HT';
        case 'Glu & 5HT'
            ntLabels{i} = 'Glu5HT';
        case 'Glu & TyramineGlu & Tyramine'
            ntLabels{i} = 'GluTyramineGluTyramine';
        case 'Glu & Tyramine'
            ntLabels{i} = 'GluTyramine';
        case 'Unknown (orphan)'
            ntLabels{i} = 'Unknown';
        end
        mask.(sprintf('%s%s',ntLabels{i},ntLabels{i})) = isType*isType';
    end

case 'AnatomyPairs'
    % Masks for all connections by anatomy type (directed)
    isHead = double(C.RegionM(:,strcmp(C.neuronAnatomyNames,'head neuron')));
    isTail = double(C.RegionM(:,strcmp(C.neuronAnatomyNames,'tail neuron')));
    isBody = double(~isHead & ~isTail);
    mask.HeadHead = isHead*isHead';
    mask.HeadTail = isHead*isTail';
    mask.TailHead = isTail*isHead';
    mask.TailTail = isTail*isTail';
    mask.BodyHead = isBody*isHead';
    mask.BodyBody = isBody*isBody';
    mask.BodyTail = isBody*isTail';
    mask.HeadBody = isHead*isBody';
    mask.TailBody = isTail*isBody';

case 'AnatomyPairsSymmetric'
    isHead = double(C.RegionM(:,strcmp(C.neuronAnatomyNames,'head neuron')));
    isTail = double(C.RegionM(:,strcmp(C.neuronAnatomyNames,'tail neuron')));
    isBody = double(~isHead & ~isTail);
    mask.HeadHead = isHead*isHead';
    mask.TailTail = isTail*isTail';
    mask.BodyBody = isBody*isBody';
    mask.HeadTail = isHead*isTail' | isTail*isHead';
    mask.BodyHead = isBody*isHead' | isHead*isBody';
    mask.BodyTail = isBody*isTail' | isTail*isBody';

case 'TypePairs'
    % Masks for all connections by type (e.g., sensory-sensory, motor-motor, etc.)
    isInter = C.RegionM(:,strcmp(C.neuronAnatomyNames,'interneuron'));
    isSensory = C.RegionM(:,strcmp(C.neuronAnatomyNames,'sensory neuron'));
    isMotor = C.RegionM(:,strcmp(C.neuronAnatomyNames,'motor neuron'));
    mask.InterInter = isInter*isInter';
    mask.InterSensory = isInter*isSensory';
    mask.SensoryInter = isSensory*isInter';
    mask.SensorySensory = isSensory*isSensory';
    mask.MotorInter = isMotor*isInter';
    mask.MotorMotor = isMotor*isMotor';
    mask.MotorSensory = isMotor*isSensory';
    mask.InterMotor = isInter*isMotor';
    mask.SensoryMotor = isSensory*isMotor';

case 'TypePairsSymmetric'
    % Masks for all connections by type (e.g., sensory-sensory, motor-motor, etc.)
    isInter = C.RegionM(:,strcmp(C.neuronAnatomyNames,'interneuron'));
    isSensory = C.RegionM(:,strcmp(C.neuronAnatomyNames,'sensory neuron'));
    isMotor = C.RegionM(:,strcmp(C.neuronAnatomyNames,'motor neuron'));
    mask.InterInter = isInter*isInter';
    mask.SensorySensory = isSensory*isSensory';
    mask.MotorMotor = isMotor*isMotor';
    mask.InterSensory = isInter*isSensory' | isSensory*isInter';
    mask.MotorInter = isMotor*isInter' | isInter*isMotor';
    mask.MotorSensory = isMotor*isSensory' | isSensory*isMotor';

case 'hubType'
    % Get hubs:
    D = GiveMeDefault();
    kHub = D.kHub;
    [~,~,deg] = degrees_dir(linkedAdj);
    isHub = double(deg > kHub)';
    maskisLR = GiveMeLRMask(C);
    ELadj = GiveMeAdj(C,'zeroBinary', 'el');
    ELadj = triu((ELadj|ELadj'),1);
    CHadj = GiveMeAdj(C,'zeroBinary', 'ch');
    CHadj = triu((CHadj|CHadj'),1);

    % Mark conns:
    connAdj = (linkedAdj|linkedAdj'); % symmetrize
    connAdj(tril(true(279))) = 0; % ignore lower diagonal

    % Link type
    mask = struct();
    mask.rich = connAdj & isHub*isHub';
    mask.feedin = connAdj & ~isHub*isHub';
    mask.feedout = connAdj & isHub*~isHub';
    mask.richfeed = mask.rich | mask.feedin | mask.feedout;
    mask.richfeedin = mask.rich | mask.feedin;
    mask.richfeedout = mask.rich | mask.feedout;

    % electrical and chemical
    mask.electrical = ELadj; 
    mask.chemical = CHadj; 
    mask.electricalOnly = ELadj & ~CHadj;  % pairs of neurons that are connected only with electrical synapses
    mask.chemicalOnly = CHadj & ~ELadj;   % pairs of neurons that are connected only with chemical synapses
    mask.chelboth = ELadj & CHadj;  % pairs of neurons that are connected with both chemical and electrical synapses
    mask.none = ~ELadj & ~CHadj; % pairs of neurons that are not connected neither with chemical, nor with electrical synapses
    mask.either = ELadj|CHadj;
    mask.CHandELonly = mask.electricalOnly|mask.chemicalOnly;

case 'unconnectedEither'
    % Ben Fulcher, 2015-02-13
    mask.special = (~linkedAdj | ~linkedAdj'); % either i->j or j->i are unconnected
    mask.special(logical(eye(size(mask.special)))) = 0; % (except self-connections)
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)

case 'unconnectedBoth'
    % Ben Fulcher, 2015-02-13
    % **both i->j and j->i
    mask.special = (~linkedAdj & ~linkedAdj'); % either i->j or j->i are unconnected
    mask.special(logical(eye(size(mask.special)))) = 0; % (except self-connections)
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)
    mask.notSpecial = (linkedAdj & linkedAdj');
    mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0; % (include just upper triangle)

    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(eye(size(mask.base)))) = 0; % (except self-connections)

case 'connectedUni'
    % Ben Fulcher, 2015-03-05
    % upper triangular unidirection connection (1-special) versus reciprocal non-connection (1-nonspecial)
    mask.special = (linkedAdj & ~linkedAdj') | (~linkedAdj & linkedAdj'); % i->j XOR j->i
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)

    mask.notSpecial = (~linkedAdj & ~linkedAdj'); % both i->j AND j->i are unconnected
    mask.notSpecial(logical(eye(size(mask.notSpecial)))) = 0; % (except self-connections)
    mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0; % (include just upper triangle)

    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(tril(ones(size(mask.base))))) = 0; % (include just upper triangle)

case 'reciprocalConnected'
    % Ben Fulcher, 2015-03-05
    % Reciprocally connected versus all others

    mask.special = (linkedAdj & linkedAdj'); % both i->j AND j->i
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)

    mask.notSpecial = ~mask.special; % not reciprocal
    mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0; % (include just upper triangle)

    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(tril(ones(size(mask.base))))) = 0; % (include just upper triangle)

case 'connectedBoth'
    % Ben Fulcher, 2015-02-13
    % upper triangular reciprocal connection (1-special) versus reciprocal non-connection (1-nonspecial)

    mask.special = (linkedAdj & linkedAdj'); % both i->j AND j->i
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)

    mask.notSpecial = (~linkedAdj & ~linkedAdj'); % both i->j AND j->i are unconnected
    mask.notSpecial(logical(eye(size(mask.notSpecial)))) = 0; % (except self-connections)
    mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0; % (include just upper triangle)

    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(tril(ones(size(mask.base))))) = 0; % (include just upper triangle)

case 'connected'
    if ~shutUp
        fprintf(1,'Differences between linked and unlinked\n');
    end
    mask.special = linkedAdj; % particular links of interest: connected
    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(eye(size(mask.base)))) = 0; % (except self-connections)

    % Not special are elements of the baseline that are not special:
    mask.notSpecial = mask.base & ~mask.special; % connections that do not exist (excluding self-connections)

case 'connectedEither'
    % Special: either connection exists:
    mask.special = (linkedAdj | linkedAdj'); % either i->j or j->i are connected
    mask.special(logical(tril(ones(size(mask.special))))) = 0; % (include just upper triangle)

    % Both unconnected:
    mask.notSpecial = (~linkedAdj & ~linkedAdj'); % both i->j AND j->i are unconnected
    mask.notSpecial(logical(eye(size(mask.notSpecial)))) = 0; % (except self-connections)
    mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0; % (include just upper triangle)

    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(tril(ones(size(mask.base))))) = 0; % (include just upper triangle)

case 'unconnected'
    % Ben Fulcher, 2015-01-15

    mask.notSpecial = linkedAdj; % connected pairs are not special
    mask.base = logical(ones(size(linkedAdj))); % baseline: everything
    mask.base(logical(eye(size(mask.base)))) = 0; % (except self-connections)

    % Not special are elements of the baseline that are not special:
    mask.special = mask.base & ~mask.notSpecial; % connections that do not exist (excluding self-connections)
% case 'ChemicalElectrical'
%     ELadj = GiveMeAdj(C,'zeroBinary', 'el');
%     CHadj = GiveMeAdj(C,'zeroBinary', 'ch');
%     mask.electricalOnly = ELadj & ~CHadj; % pairs of neurons that are connected only with electrical synapses
%     mask.chemicalOnly = CHadj & ~ELadj;  % pairs of neurons that are connected only with chemical synapses
%     mask.chelboth = ELadj & CHadj; % pairs of neurons that are connected with both chemical and electrical synapses
%     mask.none = ~ELadj & ~CHadj; % pairs of neurons that are not connected neither with chemical, nor with electrical synapses
%     mask.either = ELadj|CHadj;
case 'richfeeder'
    % ------------------------------------------------------------------------------
    % Check that there's nothing between rich and feeder:
    % Ben Fulcher, 2015-02-03
    % ------------------------------------------------------------------------------

    % Get hubs:
    kRich = extraParam{2};
    switch networkType
        case {'bu','wu'}
            nodeData = sum(linkedAdj);
        case {'bd','wd'}
            nodeData = sum(linkedAdj) + sum(linkedAdj');
    end
    isHub = (nodeData > kRich);

    % special -> rich links
    mask.special = logical(zeros(size(linkedAdj)));
    mask.special(isHub,isHub) = linkedAdj(isHub,isHub);

    % not special -> feeder links:
    mask.notSpecial = logical(zeros(size(linkedAdj)));
    mask.notSpecial(isHub,~isHub) = linkedAdj(isHub,~isHub);
    mask.notSpecial(~isHub,isHub) = linkedAdj(~isHub,isHub);

    % base --> all links (shouldn't be used)
    mask.base = linkedAdj;

case {'richlocal','feederlocal','richlocalEither','feederlocalEither'}
    % ------------------------------------------------------------------------------
    % Just rich versus local links, or feeder versus local links
    % ------------------------------------------------------------------------------
    % Ben Fulcher, 2015-01-27

    if strcmp(whatLinks,'richlocalEither') || strcmp(whatLinks,'feederlocalEither')
        % Either link exits:
        linkedAdj = (linkedAdj | linkedAdj');
    end

    % Get hubs:
    kRich = extraParam{2};
    switch networkType
        case {'bu','wu'}
            nodeData = sum(linkedAdj);
        case {'bd','wd'}
            nodeData = sum(linkedAdj) + sum(linkedAdj');
    end
    isHub = (nodeData > kRich);

    % Assign data:
    % 1. Special:
    mask.special = logical(zeros(size(linkedAdj)));
    switch whatLinks
    case {'richlocal','richlocalEither'}
        % special --> rich links
        mask.special(isHub,isHub) = linkedAdj(isHub,isHub);
    case {'feederlocal','feederlocalEither'}
        % special --> feeder links
        mask.special(isHub,~isHub) = linkedAdj(isHub,~isHub);
        mask.special(~isHub,isHub) = linkedAdj(~isHub,isHub);
    end

    % 2. Not special:
    % notSpecial --> local links
    mask.notSpecial = logical(zeros(size(linkedAdj)));
    mask.notSpecial(~isHub,~isHub) = linkedAdj(~isHub,~isHub);

    % base --> all links (not used)
    mask.base = linkedAdj;

    if strcmp(whatLinks,'richlocalEither') || strcmp(whatLinks,'feederlocalEither')
        % Just include upper triangle for either:
        mask.special(logical(tril(ones(size(mask.special))))) = 0;
        mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0;
        mask.base(logical(tril(ones(size(mask.base))))) = 0;
    end

case {'feedin','feedout','feeder','rich','richfeed','local','richfeedEither'}
    kRich = extraParam{2};

    if ~shutUp
        fprintf(1,'Differences between %s links compared to all other connected links\n',whatLinks);
    end
    switch networkType
        case {'bu','wu'}
            nodeData = sum(linkedAdj);
        case {'bd','wd'}
            nodeData = sum(linkedAdj) + sum(linkedAdj');
    end
    isRich = (nodeData > kRich);

    % If either:
    if strcmp(whatLinks,'richfeedEither')
        % Either link exits:
        linkedAdj = (linkedAdj | linkedAdj');
    end


    % Particular links of interest:
    mask.special = logical(zeros(size(linkedAdj)));
    switch whatLinks
    case 'rich'
        mask.special(isRich,isRich) = linkedAdj(isRich,isRich);
    case 'feedin'
        mask.special(~isRich,isRich) = linkedAdj(~isRich,isRich);
    case 'feedout'
        mask.special(isRich,~isRich) = linkedAdj(isRich,~isRich);
    case 'feeder'
        % Feed-in or feed-out (together):
        mask.special(isRich,~isRich) = linkedAdj(isRich,~isRich);
        mask.special(~isRich,isRich) = linkedAdj(~isRich,isRich);
    case {'richfeed','richfeedEither'}
        % Ben Fulcher, 2014-12-20
        % rich links:
        mask.special(isRich,isRich) = linkedAdj(isRich,isRich);
        % feed-in links:
        mask.special(~isRich,isRich) = linkedAdj(~isRich,isRich);
        % feed-out links:
        mask.special(isRich,~isRich) = linkedAdj(isRich,~isRich);
    case 'local'
        % Ben Fulcher, 2015-01-13
        % local links:
        mask.special(~isRich,~isRich) = linkedAdj(~isRich,~isRich);
    end

    % Baseline (all connected links)
    mask.base = linkedAdj;

    % Not special are elements of the baseline that are not special:
    mask.notSpecial = mask.base & ~mask.special;

    if strcmp(whatLinks,'richfeedEither')
        % Just include upper triangle for either:
        mask.special(logical(tril(ones(size(mask.special))))) = 0;
        mask.notSpecial(logical(tril(ones(size(mask.notSpecial))))) = 0;
        mask.base(logical(tril(ones(size(mask.base))))) = 0;
    end
end

end
