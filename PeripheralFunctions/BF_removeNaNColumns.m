% ------------------------------------------------------------------------------
% Removes columns with NaNs in them
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-09-23 Added parameter that allows some proportion of bad values
% Ben Fulcher, 2014-07-08
% ------------------------------------------------------------------------------

function [A,KeepCol] = BF_removeNaNColumns(A,propBadTol)

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2
    % Can tolerate a proportion propBad of bad values
    % Set to 0 by default (no tolerance for any bad values)
    propBadTol = 0;
end


% ------------------------------------------------------------------------------
if propBadTol==0
    KeepCol = arrayfun(@(x)~any(isnan(A(:,x))),1:size(A,2));
else
    % Compute the proportion of NaN values in each column:
    propBad = arrayfun(@(x)sum(isnan(A(:,x))),1:size(A,2))/size(A,1);
    KeepCol = (propBad > propBadTol);
end

A = A(:,KeepCol);

end