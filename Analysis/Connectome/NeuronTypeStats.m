%-------------------------------------------------------------------------------
% Simple script to output stats on neuron types (head/body/tail)
%-------------------------------------------------------------------------------
anatomyTypes = LabelNeuronType(C,'anatomy');
anatLabels = categories(anatomyTypes);
summary(anatomyTypes)

%-------------------------------------------------------------------------------
% Stats on types within anatomical regions
%-------------------------------------------------------------------------------
systemTypes = LabelNeuronType(C);
systemTypeLabels = {'interneuron','motor','sensory','multi'};
numSystemTypes = length(systemTypeLabels);
for i = 1:3
    props = arrayfun(@(x)mean(systemTypes(anatomyTypes==anatLabels{i})==systemTypeLabels{x}),1:4);
    fprintf(1,'\n%s//\n',anatLabels{i});
    for j = 1:4
        fprintf(1,'%s: %.2f%%\n',systemTypeLabels{j},props(j)*100);
    end
end

%-------------------------------------------------------------------------------
% Stats on distances:
%-------------------------------------------------------------------------------
distMat = GiveMeDist(C);
for k = 1:3
    for j = 1:3
        isTypeA = anatomyTypes==anatLabels{k};
        isTypeB = anatomyTypes==anatLabels{j};
        distances = distMat(isTypeA,isTypeB);
        maxType = max(distances(:));
        fprintf(1,'%s-%s neurons separated by at most %.2f mm\n',anatLabels{k},...
                            anatLabels{j},maxType);
    end
end
