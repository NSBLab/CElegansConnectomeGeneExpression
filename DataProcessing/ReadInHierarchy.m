function [parentName,childName,parentID,childID] = ReadInHierarchy()
% Read in hierarchy (from RetrieveHierarchy.py) and parse into matlab format

%-------------------------------------------------------------------------------
% Read in:
fid = fopen('hierarchy.csv');
hierarchy = textscan(fid, '%s%s%s%s','Delimiter','|');
fclose(fid);

%-------------------------------------------------------------------------------
% AA: 
% Changes names incerting 0 for some neurons (if name contains 3 characters and not all of them are upper letters). 
% E.g. original name in hierarchy is AS1, but listed name everywhere is AS01. This part inserts 0 where necessary. 

tf = cell(size(hierarchy{1,1},1),3);

for k=[1 3]
    tf(:,k) = isstrprop(hierarchy{:,k}, 'upper');
    for j=1:size(hierarchy{1,1},1)
        Element = hierarchy{1,k}{j};
        if length(Element)==3
            if sum(tf{j,k})<3
                hierarchy{1,k}{j} = sprintf('%s%s0%s',Element(1),Element(2),Element(3));
            end
        end
        
    end
end
%-------------------------------------------------------------------------------
% Convert to IDs:
f_getId = @(x)str2num(x(6:end));
parentName = hierarchy{1};
parentID = cellfun(f_getId,hierarchy{2});
childName = hierarchy{3};
childID = cellfun(f_getId,hierarchy{4});


%-------------------------------------------------------------------------------
% Generate translation table (parent,child):
hierarchyMap = [parentID,childID];

%-------------------------------------------------------------------------------
% Text output to user:
numRows = length(parentID);
fprintf(1,'Read data in from %s and parsed to %u hierarchical relationships\n',...
                        'hierarchy.csv',numRows);

end
