% script to check if main effect holds up with distance correction
% do correction
coExp_corr = CoExpressionDistanceCorrect(C,G); 
G.Corr.Pearson_noLR = coExp_corr;  % assign value to default coexpression measure

% run RichClub function for curves
[pConUncon, statsConUncon]=RichClub(C,G)


