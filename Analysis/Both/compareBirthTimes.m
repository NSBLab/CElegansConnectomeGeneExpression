function [coexpEarly,ph,statsh, numRich, numNONRich]=compareBirthTimes(C,G, relationship, kHub, early, connected)
% ------------------------------------------------------------------------------
% Function calculates coexpression between neurons that are born realy and
% are hubs and between those thet are born early and are not hubs.
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% G (gene data structure)
% relationship = 'born early' (neurons that are born early) or 'born40
% together' - neurons that are born together (no matter early or late -
% just together).
% khub = hub threshold (41 bu default)
% early = threshold for defining birth time for early neurons(1000min by default) or for beuront that are born together.
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute 2017-04-10
% ------------------------------------------------------------------------------%
if nargin < 3
    relationship = 'born early';
    fprintf('Do born early by DEFAULT\n')
end

if nargin < 4
    D = GiveMeDefault();
    kHub = D.kHub;
    fprintf('Hub threshold k>%u by DEFAULT\n',kHub)
end

if nargin < 5 && strcmp(relationship, 'born early')
    early = 1000;
    fprintf('Early birth time threshold is 1000min by DEFAULT\n');
elseif nargin < 5 && strcmp(relationship, 'born together')
    early = 600;
    fprintf('Born together threshold is 600min by DEFAULT\n');
end

Adj = GiveMeAdj(C,'zeroBinary');
%Adj = ones(279,279); 
[~,~,deg] = degrees_dir(Adj);
if ~connected
    Adj = ones(279);
end
isHub = (deg > kHub);
Adj=triu(Adj|Adj'); 
maskisLR = GiveMeLRMask(C); 

coexpression = GiveMeCoexpression(G,[],false); 
coexpression = coexpression.*~maskisLR; 

switch relationship
    case 'born together'
        BirthTimeDif = C.BirthTimeDiff_noLR;
        maskEarly = BirthTimeDif < early;
        maskLate = BirthTimeDif >= early;
    case 'born early'
        BirthTime = C.BirthTime;
        isEarly = BirthTime < early;
        maskEarly = isEarly.*isEarly';
        maskLate = ~isEarly.*~isEarly';
end


% take coexpression for neurons that are born early and are hubs
maskHubs = isHub.*isHub';
maskEH = (maskEarly&maskHubs&Adj)+0;
%maskEH = maskEH+NaNmatrix; 
%maskEH(isnan(maskEH))=0;
maskNONHubs = ~maskHubs +0; 
%maskEH= triu(maskEH,1); %(maskEH==0)=nan;

%maskL = (maskLate&Adj)+0;
%maskL(maskL==0)=nan;

coexpEarlyHubs = coexpression.*maskEH;
coexpEarlyHubs = nonzeros(coexpEarlyHubs); %(~isnan(coexpEarlyHubs)));

%coexpLate = coexpression.*maskL;
%coexpLate=(coexpLate(~isnan(coexpLate)));
% take coexpression between neurons that are born early and are not hubs.

% maskNONHubs2 = isHub.*~isHub';
% maskNONHubs3 = ~isHub.*isHub';
%maskNONHubs = maskNONHubs1|maskNONHubs2|maskNONHubs3;
maskENH = (maskEarly&maskNONHubs&Adj&~maskisLR)+0;
maskENH = triu(maskENH,1); 
%maskENH(maskENH==0)=nan;
%maskENH = maskENH+NaNmatrix; 
coexpEarlyNONHubs = coexpression.*maskENH;
coexpEarlyNONHubs = nonzeros(coexpEarlyNONHubs); %(coexpEarlyNONHubs(~isnan(coexpEarlyNONHubs)));

coexpEarly{1} = coexpEarlyHubs; numRich = sum(logical(coexpEarlyHubs(:)));
coexpEarly{2} = coexpEarlyNONHubs; numNONRich = sum(logical(coexpEarlyNONHubs(:)));


[ph,~,statsh] = ranksum(coexpEarlyHubs,coexpEarlyNONHubs);
%[pl,~,statsl] = ranksum(coexpEarlyNONHubs,coexpLate);

end
