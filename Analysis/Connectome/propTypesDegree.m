function propTypesDegree(C)
% ------------------------------------------------------------------------------
% Function plots the proportion of interneurons/motor/sensory/polymodal
% neurons as a function of degree.
%-------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure)
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute 2017-03-20
% ------------------------------------------------------------------------------

types = [10 13 17 21];

% chemical+electrical network
Adj = GiveMeAdj(C,'zeroBinary');

% Compute degree:
[~,~,deg] = degrees_dir(Adj);
kmax = max(deg);
kmin = min(deg);
degsort = (kmin: kmax);

% calculate proportion
prop = zeros(length(types), length(degsort));
for i = 1:length(degsort)
    ind = (deg > degsort(i));

    type = C.RegionM(ind==1,types);
    for j=1:size(type,2)
        prop(j,i) = sum(type(:,j))/sum(type(:));
    end

end

% plot proportion
f = figure('color','w');
subplot(5,1,1)
kRange = degsort;
N = arrayfun(@(x)sum(deg==x),kRange);
bar(kRange,N,'EdgeColor','k','FaceColor','k')
xlim([min(deg)-0.5,max(deg)+0.5]);
ylabel('Frequency')
subplot(5,1,2:5)
% Give colors:
colors = BF_getcmap('spectral',length(types),1);
kRange = min(deg):max(deg);
for i = 1:length(kRange)
    ind = (deg > degsort(i));
    categoriesHere = C.RegionM(ind==1, types);
    proportions = prop'; %countcats(categoriesHere)/length(categoriesHere);
    for j = 1:size(proportions,2)
        if proportions(i,j) > 0
            rectangle('Position',[kRange(i),sum(proportions(i,1:j-1)),1,proportions(i,j)], ...
                'FaceColor',colors{j},'EdgeColor',colors{j})
        end
    end
end
Labels = C.neuronAnatomyNames(types);
where = zeros(length(types),1);
for w=1:length(types)
    if w==1
        where(w) = proportions(1,w)/2;
    else
        where(w) = sum(proportions(1,1:w-1))+proportions(1,w)/2;
    end
end

prop1 = proportions(1,:);
for j = 1:length(prop1)
    text(kRange(1),where(j),Labels{j})
end
end
