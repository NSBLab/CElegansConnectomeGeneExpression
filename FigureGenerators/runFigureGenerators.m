% ------------------------------------------------------------------------------
% Aurina Arnatkeviciute, 2017-03-22
% ------------------------------------------------------------------------------
% script to run functions that generate plots

%% normalised RC coefficient
networkRC(C)

%% plot coexpression between connected and unconnected links
[p1, stats1] = coexpConnUncon(C, G);

%% plot coexpression (or lineage/connection_distance) for different types of links using mean to sumarise a measure at each degree threshold
RichClub(C,G,'coexpression'); %ylim([0 0.7]);
RichClub(C,G,'lineage distance'); ylim([14 20]);
RichClub(C,G,'connection distance'); ylim([0 1]);

%% plot coexpression (or lineage/connection_distance) for different types of links using mean to sumarise a measure at each degree threshold
RichClubMedian(C,G,'coexpression');  ylim([0.9 1]);
% RichClubMedian(C,G,'lineage distance'); % doesn't make sense;
RichClubMedian(C,G,'connection distance'); ylim([0 1.2]);

%% plot coexpression as a function of degree for interneurons, sensory neurons and motor neurons.
[data, numberHere] = degreeDependence(C,G);

%% plot coexpression for hub interneurons and non-hub interneurons
[p2, stats2] = plotInterneuronsDegree(C,G);

%% modules
gam = 1;
for gamma=gam
    [pwithinbetween, z, pwithinbetweenNOR, statswithinbetweenNOR, M2, BTWstats]=communityCoexpression(C,G, gamma);
end
