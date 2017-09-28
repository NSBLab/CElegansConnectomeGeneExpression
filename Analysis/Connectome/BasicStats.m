%-------------------------------------------------------------------------------
% Get basic stats on connections
%-------------------------------------------------------------------------------
Adj_w = GiveMeAdj(C,'zeroWeighted');
Adj_b = GiveMeAdj(C,'zeroBinary');
Adj_b_sym = double(triu(Adj_b|Adj_b',+1));

fprintf(1,'%u weighted directed connections\n',sum(Adj_w(:)));
fprintf(1,'%u binary directed connections\n',sum(Adj_b(:)));
fprintf(1,'%u binary pairwise connections\n',sum(Adj_b_sym(:)));

N = length(Adj_b_sym);
fprintf(1,'%.2f%% neuron pairs are connected\n',100*sum(Adj_b_sym(:))/(N*(N-1)/2));
