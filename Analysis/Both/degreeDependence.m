%%  [data, numberHere] = degreeDependence(C,G, measureTouse, doStrength, conType, chooseCategory,  includingCategory, chooseRigth, useKinKout, doNonOverlap, categoryNoOverlap)
% ------------------------------------------------------------------------------
% function plots coexpression for links involving a specific category
% (motor/sensory/interneurons) as a function of degree (or strength).
% ------------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------------
% C (connectivity structure),
% G (gene data structure);
% coexpMeasure - choose coexpression measure as defined in G.Corr. default Pearson_noLR
% measureTouse - mean or median to summarise coexpression at each threshold
% conType - chemical (ch), all(all), electrical(el)
% chooseRight - 0/1; on top of the category chosen, it will filter only
% right side neurons of that category;
% chooseCategory = [] - set of numbers from C.C.neuronAnatomyNames (can be from 1 to 10)
% includingCat = 0/1/2/3
% (0 - exclude connections involving neurons in subcategory;
%  1 - include connections involving neurons in subcategory,
%  2 - include only connections within subcategory (subcategory neurons
%  connecting to other subcategory neurons)
%  3 - exclude connections involving neurons within subcategory (e.g. include all connections except interneuron-interneuron connections)
% ------------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------------
% data - degree, coexpression and number fo links for each type of neurons;
% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2017-03-20
% ------------------------------------------------------------------------------
function [data, numberHere] = degreeDependence(C,G, coexpMeasure, measureTouse, doStrength, conType, chooseCategory,  includingCategory, chooseRigth, useKinKout, doNonOverlap)

if nargin < 3
    correlation = GiveMeCoexpression(G,[],true);
end

if nargin < 4
    measureTouse = 'mean';
end

if nargin < 5
    doStrength = false;
end
if nargin < 6
    conType = 'all';
end

if nargin < 7
    chooseCategory = [10 13 17]; %interneuron, motor, sensory;
end

if nargin < 8
    includingCategory = 1;
end

if nargin < 9
    chooseRigth = false;
end

if nargin < 10
    useKinKout = false;
end

if nargin < 11
    doNonOverlap = false;
end


Adj = GiveMeAdj(C,'zeroWeighted',conType);
correlation = GiveMeCoexpression(G,coexpMeasure,true); 

if doStrength
    [kin,kout,deg] = strengths_dir(Adj);
else
    [kin,kout,deg] = degrees_dir(Adj);
end


if useKinKout == 1
    deg = kout;
end
Adj = logical(Adj);
Adj=logical(Adj+Adj')+0;
nans = tril(nan(279,279)); % make nan matrix and use that to exclude the lower triangle.
Adj = Adj+nans;

coexp = correlation.*Adj;
coexp(coexp==0) = NaN;


data = cell(length(chooseCategory),2);
figure('color', 'w');
set(gcf,'Position',[1000,200,1000,500])
n=ones(1,279);
for k=1:length(chooseCategory)
    l = chooseCategory(k);
    subCat = C.RegionM(:,l)';
    %     if useOppositeCat == 1
    %         subCat = ~subCat;
    %     end
    if doNonOverlap == 1

        %subCat(subCat~=3)=0; - interneurons, but not in the head
        %subCat(subCat~=1)=0; - head neurons, but not interneurons
        subCat(subCat~=1)=0;
        subCat = logical(subCat);
    end
    if chooseRigth == 1

        isRight = extractfield(C.RegionStruct, 'Right');
        subCat = subCat+isRight;
        thr = 2;

    else
        thr = 1;
    end

    [uniqueDeg, expHere,numSubcat] = coexpressionDeg(subCat, coexp, includingCategory);
    lineStyle = '-'; markerStyle = 'o';
    data{k,1} = uniqueDeg;
    data{k,2} = expHere;
    data{k,3} = numSubcat;
    uniqueDegrees = unique(horzcat(data{:,1}));
    colorlist = GiveMeColors('InterneuronMotorSensoryMulti'); %('[1 .46 .22; .32 .28 .31; 0 .5 .5];%{'r', 'b', 'c', 'm', 'g', 'k', 'm'};%'color', CM(k,:),
    sp2=subplot(2,5,6:8); plot(data{k,1}, data{k,2}, lineStyle, 'color',colorlist(k,:),'LineWidth',2.5); hold on;
    plot(data{k,1}, data{k,2},markerStyle,'MarkerEdgeColor',colorlist(k,:),...
            'MarkerFaceColor',brighten(colorlist(k,:),+0.5),'LineWidth',1,'MarkerSize',7);

    legendInfo{k} = [C.neuronAnatomyNames{l}];box off;

end

pos=get(sp2,'Position');
    set(sp2,'Position',[pos(1)*1, pos(2)*2.5, pos(3)*1, pos(4)*0.8]); % [left bottom width height]
%legend(legendInfo);
% if includingCategory == 1
%     title('Connections involving category');
% elseif includingCategory == 0
%     title('Connections not involving category');
% elseif includingCategory == 2
%     title('Connections only between category');
% elseif includingCategory == 3
%     title('Connections only between non category');

%end
if doStrength == 1
    xlabel('strength>=s', 'FontSize', 12);
else
    xlabel('Degree, k', 'FontSize', 12);
end
switch measureTouse
    case 'mean'
        ylabel('Gene coexpression', 'FontSize', 12);
    case 'median'
        ylabel('Gene coexpression', 'FontSize', 12);
end

uniqueDegrees = unique(horzcat(data{:,1}));
numberHere = zeros(length(uniqueDegrees), length(chooseCategory));
for p=1:length(uniqueDegrees)
    for u=1:length(chooseCategory)
        dcat = uniqueDegrees(p);
        indHere= find(data{u,1}==dcat);
        if isempty(indHere)
            numberHere(p,u) = 0;
        else
            numberHere(p,u) = data{u,3}(indHere);
        end
    end


end

sp1=subplot(2,5,1:3); H = bar(uniqueDegrees,numberHere, 'stacked', 'EdgeColor', 'w'); xlim([1,max(uniqueDegrees)]); ylim([0,3500]);
set(gca,'Xtick', []);
pos=get(sp1,'Position');
    set(sp1,'Position',[pos(1)*1, pos(2)*1, pos(3)*1, pos(4)*0.5]); % [left bottom width height]
%title('Number of links for each category');

if doStrength == 1
    ylabel('Number of links', 'FontSize', 12);%xlabel('strength>=s');
else
    ylabel('Number of links', 'FontSize', 12);%xlabel('degree>=k');
end
for v=1:length(chooseCategory)
    set(H(v),'facecolor',colorlist(v,:));
end
box off; hold on;
%AX=legend(H, {'a','b','c','d','e','f'}, 'Location','Best','FontSize',8);
    function [uniqueDeg, expHere, numSubcat] = coexpressionDeg(subCat, coexp,includingCategory)
        if includingCategory==1

            nc1 = sum(subCat); valsc1=linspace(1,nc1,nc1); c1rank = zeros(279,1); c1rank(subCat==1)=valsc1; valsc12 = linspace(nc1+1,279,279-nc1); c1rank(subCat==0)=valsc12;
            [~, valInd] = sort(c1rank);

            categoryDeg = deg(valInd);
            category = coexp(valInd, valInd);
            category((nc1+1):279, (nc1+1):279)=NaN;
            %             calcCatDeg = categoryDeg(1:nc1);

        elseif includingCategory==0

            nc1 = sum(~subCat); valsc1=linspace(1,nc1,nc1); c1rank = zeros(279,1); c1rank(subCat==0)=valsc1; valsc12 = linspace(nc1+1,279,279-nc1); c1rank(subCat==1)=valsc12;
            [~, valInd] = sort(c1rank);

            categoryDeg = deg(valInd);
            category = coexp(valInd, valInd);
            category((nc1+1):279, (nc1+1):279)=NaN;
            %             calcCatDeg = categoryDeg(1:nc1);

        elseif includingCategory==2

            category = coexp(subCat == thr,subCat == thr);
            categoryDeg = deg(subCat==thr);

        elseif includingCategory==3

            subCat=~subCat;
            category = coexp(subCat == thr,subCat == thr);
            categoryDeg = deg(subCat==thr);

        end

        uniqueDeg = unique(categoryDeg);
        expHere = zeros(length(uniqueDeg),1);
        numSubcat = zeros(length(uniqueDeg),1);
        for j=1:length(uniqueDeg)
            degHere = uniqueDeg(j);
            indDeg = find(categoryDeg>degHere);
            %calcCatDeg = categoryDeg(
            %numSubcat(j) = length(indDeg);
            numLinks = category(indDeg,indDeg);
            numLinks(isnan(numLinks))=0;
            numSubcat(j) = sum(sum(logical(numLinks)));
            switch measureTouse
                case 'median'
            expHere(j) = nanmedian(nanmedian(category(indDeg,indDeg)));
                case 'mean'
            expHere(j) = nanmean(nanmean(category(indDeg,indDeg)));
            end

        end
        numSubcat(isnan(expHere))=NaN;
        uniqueDeg(isnan(expHere))=NaN;

    end
% plot non-hub/ hub interneurons
[p, stats, dataCell] = plotInterneuronsDegree(C,G);
hold on;
sp3 = subplot(2,5,4:5);

extraParams = struct('customSpot','');
JitteredParallelScatter(dataCell,0,1,false, extraParams);
xLabels = {'Hub interneurons', 'Non-hub interneurons'};
set(gca,'Xtick', [1 2], 'XTickLabel',xLabels, 'FontSize', 12);
set(gca,'Ytick', [0 0.2 0.4 0.6 0.8 1], 'YTickLabel',[0 0.2 0.4 0.6 0.8 1], 'FontSize', 12);
set(gca,'box','off');
ylabel('Coexpression','FontSize', 12);
pos=get(sp3,'Position');
    set(sp3,'Position',[pos(1)*1.05, pos(2)*0.5, pos(3)*1, pos(4)*1.2]); % [left bottom width height]

end
