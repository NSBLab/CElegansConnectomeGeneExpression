% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-09-18
% ------------------------------------------------------------------------------
% Quick implementation of communicability, a measure introduced by Estrada and
% Hatano (2008).
% 
%---INPUTS:
% Adj, the adjacency matrix

function G = communicability(Adj)

% Check that diagonal is zero
if ~all(diag(Adj)==0)
	fprintf(1,'Setting diagonal to zero.\n');
	Adj(logical(eye(size(Adj)))) = 0;
end

if islogical(Adj)
    Adj = double(Adj);
end

% ------------------------------------------------------------------------------
% Compute the communicability measure as the exponential of the adjacency matrix.
% ------------------------------------------------------------------------------
% This weights higher order paths between two nodes lower, decreasing as 1/x!
G = expm(Adj);

% Set diagonal to zero:
% G(diag(G)) = 0;

% ------------------------------------------------------------------------------
% Compute mean communicability, meanG:
% ------------------------------------------------------------------------------
% (keeping diagonal terms in the mix when computing the average):
% meanG = mean(G(logical(Adj)));

end