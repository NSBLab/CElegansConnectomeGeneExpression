function PlotColorMatrix(DataMatrix,RegionStruct,ColorLabelsWhere,rectThickness,MyColorMap,labelInd,isPermuted,extraParams)
% Plots a colored data matrix, with the mouse connectome regions labeled
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-03
% Modified by Taqi Ali 20-7-15
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 3
    fprintf(1,'We need colors and labels for each element of the matrix\n');
end
if nargin < 3 || isempty(ColorLabelsWhere)
    ColorLabelsWhere = 'left'; % left, bottom, both
end
if nargin < 4 || isempty(rectThickness)
    rectThickness = arrayfun(@(x)size(DataMatrix,x)/50,1:2);
end
if length(rectThickness)==1, rectThickness = ones(2,1)*rectThickness; end
if nargin < 5 || isempty(MyColorMap)
    MyColorMap = flipud(BF_getcmap('redblue',11,0));
    % MyColorMap = flipud(BF_getcmap('redyellowblue',9,0));
end
if nargin < 6 || isempty(labelInd)
    % Labels individual regions rather than major regions
    labelInd = 0;
end
if nargin < 7 || isempty(isPermuted)
    isPermuted = 0; % Assume it's in the region order
end
if nargin < 8
    extraParams = struct;
end

% ------------------------------------------------------------------------------
% Check extraParams
if ~isfield(extraParams,'plotBoundaries')
    extraParams.plotBoundaries = (size(DataMatrix)==213);
    % plot boundaries if dimensions are 213
end

% ------------------------------------------------------------------------------
% Provided just a single RegionStruct -- set the same for rows and columns
% ------------------------------------------------------------------------------
if ~iscell(RegionStruct)
    RegionStruct_in = RegionStruct;
    RegionStruct = cell(2,1);
    RegionStruct{1} = RegionStruct_in;
    RegionStruct{2} = RegionStruct_in;
end

% Flipped version because pcolor is weird:
RegionStruct_flip = RegionStruct{1}(length(RegionStruct{1}):-1:1); % first labels at the top

% ------------------------------------------------------------------------------
% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
DataMatrix = flipud(DataMatrix); % first labels at the top
pcolor([DataMatrix, zeros(size(DataMatrix,1),1); zeros(1,size(DataMatrix,2)+1)]);
shading flat
colormap(MyColorMap)
hold on

% ------------------------------------------------------------------------------
% Superimpose black rectangles over NaN values
% ------------------------------------------------------------------------------
% if any(isnan(DataMatrix(:)))
%     if ~isfield(extraParams,'NaNcolor')
%         extraParams.NaNcolor = 'k';
%     end
%     Gray = ones(3,1)*0.7;
%     [theNaNs_i,theNaNs_j] = find(isnan(DataMatrix));
%     for i = 1:length(theNaNs_i)
%         rectangle('Position',[theNaNs_j(i),theNaNs_i(i),1,1],'FaceColor',extraParams.NaNcolor, ...
%                         'EdgeColor',extraParams.NaNcolor)
%     end
% end
if ~isfield(extraParams,'backgroundColor')
    extraParams.backgroundColor = 'k';
end
set(gca,'color',extraParams.backgroundColor);

% ------------------------------------------------------------------------------
% Y-axis labels in the middle of each contiguous region of major region labels
% ------------------------------------------------------------------------------
% Get major region labels, and their positions (for text labels, also squares later):
Labels = {RegionStruct_flip.acronym};
[~,ia,ib] = unique(Labels,'stable');

if labelInd || isPermuted
    % Label individual regions:
    set(gca,'YTick',1.5:length(RegionStruct{1})+1);
    % set(gca,'YTickLabel',{RegionStruct_flip.acronym});
    set(gca,'YTickLabel',arrayfun(@(x)sprintf('%s (%s)',RegionStruct_flip(x).name,...
                        RegionStruct_flip(x).acronym),1:length(RegionStruct{1}),'UniformOutput',0));

    % X-tick-labels
    set(gca,'XTick',1.5:length(RegionStruct{2})+1);
    set(gca,'XTickLabel',arrayfun(@(x)sprintf('%s',RegionStruct{2}(x).acronym),1:length(RegionStruct{2}),'UniformOutput',0));
else
    % Label major region names:
    ia_mid = [ia;length(Labels)];
    ia_mid = floor(mean([ia_mid(1:end-1),ia_mid(2:end)],2));
    set(gca,'YTick',ia_mid);
    set(gca,'YTickLabel',Labels(ia_mid));
end

% ------------------------------------------------------------------------------
% Add squares for each region in the case of square matrix; separator lines in
% the case of a rectangular matrix
% ------------------------------------------------------------------------------
if (any(extraParams.plotBoundaries) && ~isPermuted)
    LineWidth = 2;

    if size(DataMatrix,1) == size(DataMatrix,2) % a square matrix
        iap = ia;
        iap(end) = iap(end);
        N = size(DataMatrix,1);
        iadiff = diff([iap; length(Labels)+1]);
        for i = 1:length(ia)
            ColorHere = rgbconv(RegionStruct_flip(ia(i)).color_hex_triplet);
            rectangle('Position',[N-iap(i)-iadiff(i)+2,iap(i),iadiff(i),iadiff(i)], ...
                            'EdgeColor',ColorHere,'LineWidth',LineWidth,'LineStyle','-')
        end
        axis square

    else

        % Separator lines---rows
        if extraParams.plotBoundaries(1)

            for i = 1:length(ia)
                ColorHere = rgbconv(RegionStruct_flip(ia(i)).color_hex_triplet);
                plot([0.5,size(DataMatrix,2)+1.5],ia(i)*ones(2,1),'--', ...
                                'color',ColorHere,'LineWidth',LineWidth)
            end
        end

        % Separator lines---columns
        if extraParams.plotBoundaries(2)
            [~,ia_col,~] = unique({RegionStruct{2}.MajorRegionName},'stable');
            for i = 1:length(ia_col)
                ColorHere = rgbconv(RegionStruct{2}(ia_col(i)).color_hex_triplet);
                plot(ia_col(i)*ones(2,1),0.5+[0,size(DataMatrix,1)+1],'--', ...
                                'color',ColorHere,'LineWidth',LineWidth)
            end
        end
    end
end

% ------------------------------------------------------------------------------
% Add rectangles labeling major brain regions, and individual colors
% ------------------------------------------------------------------------------

% Rows:
% for j = 1:length(RegionStruct{1})
%     % My colors:
%     % rectangle('Position',[-2*rectThickness,j,rectThickness,1],'FaceColor',c{ib(j)},'EdgeColor',c{ib(j)})
%
%     % Add rectangle to color each row (perhaps also column):
%
%     if ismember(ColorLabelsWhere,{'both','left'})
%         ColorHere = rgbconv(RegionStruct_flip(j).color_hex_triplet);
%         rectangle('Position',[1-rectThickness(2),j,rectThickness(2),1], ...
%                     'FaceColor',ColorHere,'EdgeColor',ColorHere)
%     end
%     if ismember(ColorLabelsWhere,{'both','right'})
%         ColorHere = rgbconv(RegionStruct_flip(j).color_hex_triplet);
%         rectangle('Position',[size(DataMatrix,2)+1,j,rectThickness(2),1], ...
%                     'FaceColor',ColorHere,'EdgeColor',ColorHere)
%     end
% end

% Columns:
% for j = 1:length(RegionStruct{2})
%     if ismember(ColorLabelsWhere,{'both','bottom'})
%         ColorHere = rgbconv(RegionStruct{2}(j).color_hex_triplet);
%         rectangle('Position',[j,1-rectThickness(1),1,rectThickness(1)], ...
%                     'FaceColor',ColorHere,'EdgeColor',ColorHere)
%     end
%     if ismember(ColorLabelsWhere,{'both','top'})
%         ColorHere = rgbconv(RegionStruct{2}(j).color_hex_triplet);
%         rectangle('Position',[j,size(DataMatrix,1)+1,1,rectThickness(1)], ...
%                     'FaceColor',ColorHere,'EdgeColor',ColorHere)
%     end
% end

% ------------------------------------------------------------------------------
% Add separator black lines:
% ------------------------------------------------------------------------------
LineWidth = 1.2;
if ismember(ColorLabelsWhere,{'both','left'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot([1-rectThickness(2),1],ones(2,1)*ia(j),'k','LineWidth',LineWidth)
        end
    end
    % Bottom one:
    plot([1-rectThickness(2),1],ones(2,1)*1,'k','LineWidth',LineWidth)
    % Top one:
    plot([1-rectThickness(2),1],ones(2,1)*size(DataMatrix,1)+1,'k','LineWidth',LineWidth)
    % Left one:
    plot((1-rectThickness(2))*ones(2,1),[1,size(DataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Right one:
    plot(ones(2,1),[1,size(DataMatrix,1)+1],'k','LineWidth',LineWidth)
end
if ismember(ColorLabelsWhere,{'both','right'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(size(DataMatrix,2)+1+[0,rectThickness(2)],ones(2,1)*ia(j),'k','LineWidth',LineWidth)
        end
    end
    % Top one:
    plot(size(DataMatrix,2)+1+[0,rectThickness(2)],ones(2,1)*size(DataMatrix,1)+1,'k','LineWidth',LineWidth)
    % Left one:
    plot((size(DataMatrix,2)+1)*ones(2,1),[1,size(DataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Right one:
    plot((size(DataMatrix,2)+1+rectThickness(2))*ones(2,1),[1,size(DataMatrix,1)+1],'k','LineWidth',LineWidth)
    % Bottom one:
    plot(size(DataMatrix,2)+1+[0,rectThickness(2)],ones(2,1),'k','LineWidth',LineWidth)
end
if ismember(ColorLabelsWhere,{'both','bottom'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(ones(2,1)*(size(DataMatrix,1)-ia(j)+2),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
        end
    end
    % Leftmost one:
    plot(ones(2,1),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
    % Bottom one:
    plot([1,size(DataMatrix,2)+1],(1-rectThickness(1))*ones(2,1),'k','LineWidth',LineWidth)
    % Right one:
    plot(ones(2,1)*(size(DataMatrix,2)+1),[1-rectThickness(1),1],'k','LineWidth',LineWidth)
    % Top one:
    plot([1,size(DataMatrix,2)+1],ones(2,1),'k','LineWidth',LineWidth)
end
if ismember(ColorLabelsWhere,{'both','top'})
    if ~(labelInd || isPermuted)
        for j = 1:length(ia)
            plot(ones(2,1)*(size(DataMatrix,2)-ia(j)+2),size(DataMatrix,1)+1+[0,rectThickness(1)], ...
                            'k','LineWidth',LineWidth)
        end
    end
    % Leftmost one:
    plot(ones(2,1),size(DataMatrix,1)+1+[0,rectThickness(1)],'k','LineWidth',LineWidth)
    % Rightmost one:
    plot(ones(2,1)*(size(DataMatrix,2)+1),size(DataMatrix,1)+1+[0,rectThickness(1)],'k','LineWidth',LineWidth)
    % Bottom one:
    plot([1,size(DataMatrix,2)+1],(size(DataMatrix,1)+1)*ones(2,1),'k','LineWidth',LineWidth)
    % Top one:
    plot([1,size(DataMatrix,2)+1],(size(DataMatrix,1)+1+rectThickness(1))*ones(2,1),'k','LineWidth',LineWidth)
end

% ------------------------------------------------------------------------------
% Adjust axes to see the labeling:
% ------------------------------------------------------------------------------
scaleFactor = 0.08; % to see the little bit extra to capture the line thickness
% First set with the scale factor
switch ColorLabelsWhere
case 'both'
    set(gca,'XLim',[-rectThickness(2)*(scaleFactor),size(DataMatrix,2)+rectThickness(2)*scaleFactor])
    set(gca,'YLim',[-rectThickness(1)*(scaleFactor),size(DataMatrix,1)+rectThickness(1)*scaleFactor])
case 'left'
    set(gca,'XLim',[-rectThickness(2)*(scaleFactor),size(DataMatrix,2)+rectThickness(2)*scaleFactor])
end
if ismember(ColorLabelsWhere,{'both','left'})
    % get_xlim = get(gca,'xlim');
    % get_xlim(1) = -rectThickness*(1+scaleFactor);
    set(gca,'XLim',[1-rectThickness(2)*(1+scaleFactor),size(DataMatrix,2)+1]);
end
if ismember(ColorLabelsWhere,{'both','right'})
    get_xlim = get(gca,'xlim');
    get_xlim(2) = size(DataMatrix,2) + 2 + rectThickness(2)*(1+scaleFactor);
    set(gca,'XLim',get_xlim)
end
if ismember(ColorLabelsWhere,{'both','bottom'})
    get_ylim = get(gca,'ylim');
    get_ylim(1) = -rectThickness(1)*(1+scaleFactor);
    set(gca,'YLim',get_ylim)
end
if ismember(ColorLabelsWhere,{'both','top'})
    get_ylim = get(gca,'ylim');
    get_ylim(2) = size(DataMatrix,1) + 2 + rectThickness(1)*(1+scaleFactor);
    set(gca,'YLim',get_ylim)
end

% Remove tick marks:
set(gca,'TickLength',[0,0])

end
