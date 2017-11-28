%===============================================================================
% Coexpression for different types of links: use network comprised of
% chemical and electrical synapses
%===============================================================================
[S,P] = RichClub(C,G, 'all');
%===============================================================================
% Coexpression for different types of links: use network comprised only of
% chemical synapses
%===============================================================================
[S,P] = RichClub(C,G, 'synaptic');