function maskIsLR = GiveMeLRMask(C)
% Outputs a mask which is true where there are LR matches
%-------------------------------------------------------------------------------
RegionStruct = C.RegionStruct; 
leftNeuronIx = find([RegionStruct.isLeft]);
% Make mask for removing coexpression for L/R pairs of neurons
maskIsLR = false(279,279);

for i=1:length(leftNeuronIx)
    index_L = leftNeuronIx(i);
    nameBase = RegionStruct(index_L).acronym(1:end-1); % for each left neuron find right neuron, check if names match and remove coexpression value
    NameR = [nameBase, 'R'];
    index_R = find(strcmp({RegionStruct.acronym},NameR));
    if isempty(index_R)
        warning('No right-match for %s',RegionStruct(index_L).acronym)
        continue
    end
    % find matching right
    if strcmp(RegionStruct(index_L).acronym(end),'L') && ...
            strcmp(RegionStruct(index_R).acronym(end),'R') && ...
            strcmp(RegionStruct(index_L).acronym(1:end-1),RegionStruct(index_R).acronym(1:end-1))
        % NaN the index_L/index_R entry as a symmetric match
        maskIsLR(index_L,index_R) = true;
        maskIsLR(index_R,index_L) = true;
    else
        warning('---Error matching left-right neurons for %s and %s',...
                    RegionStruct(index_L).acronym,RegionStruct(index).acronym);
    end
end

% Check nothing on diagonal:
if any(diag(maskIsLR))
    fprintf(1,'Diagonal entries in mask?\n');
    maskIsLR(eye(true(size(maskIsLR)))) = false;
    % maskIsLR(1:size(maskIsLR)+1:end) = NaN;  % remove diagonal
end

end
