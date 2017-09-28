%% Taqi Ali 23-7-15
% Plots the degree distribution

%% Retrieve Degree information

% Binary Ch
Adj_ch = GiveMeAdj(C,'ch','zeroBinary');
[~,~,deg_ch] = degrees_dir(Adj_ch);

% Binary gap junctions
Adj_ej = GiveMeAdj(C,'ej','zeroBinary');
deg_ej = degrees_und(Adj_ej);

% Binary All
Adj_all = GiveMeAdj(C,'all','zeroBinary');
[~,~,deg_all_d] = degrees_dir(Adj_all);

%-------------------------------------------------------------------------------
uniquelist = C.MajorRegionLabels;

figure('color','w')
histogram(deg_ch');
title('Directed Chemical Synapse Degree Histogram')
xlabel('k = kin + kout');
ylabel('frequency');

figure('color','w')
bar(deg_ch');
title('Directed Chemical Synapse Degree Bargraph')
set(gca,'XTick',1:length(uniquelist))
set(gca,'XTickLabel',uniquelist)
ylabel('k = kin + kout');

figure('color','w')
histogram(deg_ej');
title('Undirected Electrical Synapse Degree Histogram')
xlabel('k total');
ylabel('frequency');

figure('color','w')
bar(deg_ej');
title('Undirected Electrical Synapse Degree Bargraph')
set(gca,'XTick',1:length(uniquelist))
set(gca,'XTickLabel',uniquelist)
ylabel('k total');

figure('color','w')
histogram(deg_all_d');
title('Directed Chemical + Electrical Synapse Degree Histogram')
xlabel('k = kin + kout');
ylabel('frequency');

figure('color','w')
bar(deg_all_d');
title('Directed Chemical + Electrical Synapse Degree Bargraph')
set(gca,'XTick',1:length(uniquelist))
set(gca,'XTickLabel',uniquelist)
ylabel('k total');
