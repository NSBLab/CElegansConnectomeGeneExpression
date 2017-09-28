% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-24
% ------------------------------------------------------------------------------

function PlotCorrMat(D_corr,Minus1To1)

if nargin < 2 || isempty(Minus1To1)
    Minus1To1 = 1; % assume range from -1 to 1
end

% Make a matrix
if any(size(D_corr)==1)
    D_corr = squareform(D_corr);
end

imagesc(D_corr)

axis square

% Set color limits and colormap
if Minus1To1
    caxis([-1,1])
    colormap([flipud(BF_getcmap('blues',9,0));BF_getcmap('reds',9,0)])
else % assume [0,1] (a normalized distance metric)
    caxis([0,1])
    colormap(BF_getcmap('reds',9,0))
end

% ------------------------------------------------------------------------------
% Superimpose green/yellow rectangles over NaN values
% ------------------------------------------------------------------------------
if any(isnan(D_corr(:)))
    Green = BF_getcmap('greens',3,1);
    Red = BF_getcmap('reds',3,1);
    [theNaNs_i,theNaNs_j] = find(isnan(D_corr));
    fprintf(1,['Superimposing green/red rectangles over all %u NaNs in ' ...
                                'the data matrix\n'],length(theNaNs_i));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor',Green{end}, ...
                        'EdgeColor',Red{end})
    end
end

% Add a colorbar
colorbar

% Remove ticks:
set(gca,'TickLength',[0,0])

end