% Removes NaN entries from an input distance matrix, R
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-25
% ------------------------------------------------------------------------------
function [R, keepers] = RemoveNaN_DistMat(R)

keepers = logical(ones(length(R),1));

if any(isnan(R(:)))
    AreNaN = 1;

    while AreNaN
        NumNaNs = sum(isnan(R(keepers,keepers)));
        [~,irem] = max(NumNaNs); % the index (of keepers==1) that has the most NaNs
        fkeep = find(keepers);
        keepers(fkeep(irem)) = 0; % remove this index from the full list, keepers
        R_keep_tmp = R(keepers,keepers);
        AreNaN = any(isnan(R_keep_tmp(:))); % are there still more NaNs after removing this?
    end
    R = R(keepers,keepers);
end
 
end