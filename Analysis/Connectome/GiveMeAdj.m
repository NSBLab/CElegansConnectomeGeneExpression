function theAdjMat = GiveMeAdj(C,whatAdj,whatType)
% Gives a string identifying the type of normalization to apply, then returns
% the gene data for that normalization.
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-07-17
% Modified by Taqi Ali, 23/7/15
% whatAdj = 'rawBinary','raw weighted','zero binary','zero weighted'
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(whatAdj)
    whatAdj = 'zeroWeighted';
end
if nargin < 3 || isempty(whatType)
    whatType = 'all'; %'ch','ej','all'
end

% ------------------------------------------------------------------------------
% Take the appropriate connections
% ------------------------------------------------------------------------------

switch whatType
case 'ch' % chemical synapses
    theIndex = 1;
case 'el' % electrical gap junctions
    theIndex = 2;
case 'all' % both chemical and electrical connections
    theIndex = 3;
end

% ------------------------------------------------------------------------------
% Retrieve and process the data
% ------------------------------------------------------------------------------
switch whatAdj
case 'zeroWeighted'
    % Raw weights but with zeros where either zero or NaN
    theAdjMat = C.Adj_W{theIndex};

    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;

case 'rawWeighted'
    % Raw weighted
    theAdjMat = C.Adj_W{theIndex};

case 'zeroBinary'
    % binary but with zeros where either zero or NaN
    theAdjMat = C.Adj_B{theIndex};

    % Remove diagonal entries:
    theAdjMat(logical(eye(size(theAdjMat)))) = 0;

case 'rawBinary'
    % rawBinary
    theAdjMat = C.Adj_B{theIndex};

otherwise
    error('Unknown adjacency matrix option: ''%s''',whatAdj);
end

end
