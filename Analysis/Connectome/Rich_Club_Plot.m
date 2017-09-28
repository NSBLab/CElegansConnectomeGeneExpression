%% Taqi Ali 23-7-15
% Creates rich-club coefficient plot for non-normalised data
%-------------------------------------------------------------------------------

[R_ch, Nk_ch, Ek_ch] = rich_club_bd(C.Adj_B{1});
figure
plot(1:length(R_ch),R_ch)
title('Rich-Club Coefficient for Directed Chemical Synapses');
xlabel('k (degree)');
ylabel('RCC \phi(k)');

[R_ej, Nk_ej, Ek_ej] = rich_club_bu(C.Adj_B{2});
figure
plot(1:length(R_ej),R_ej)
title('Rich-Club Coefficient for Undirected Electrical Synapses');
xlabel('k (degree)');
ylabel('RCC \phi(k)');

[R_all_d, Nk_all_d, Ek_all_d] = rich_club_bd(C.Adj{3});
figure
plot(1:length(R_all_d),R_all_d)
title('Rich-Club Coefficient for Directed Chemical + Electrical Synapses');
xlabel('k (degree)');
ylabel('RCC \phi(k)');
