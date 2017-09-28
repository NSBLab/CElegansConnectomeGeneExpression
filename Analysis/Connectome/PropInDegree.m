% Investigating the proportion of links that are incoming as a function of total degree
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-09-30
% ------------------------------------------------------------------------------

pThreshold = 1.1;

Adj = GiveMeAdj(C,'binary','ipsi',0,pThreshold);

inDegree = sum(Adj,1);
totDegree = sum(Adj,1) + sum(Adj,2)';

nodeData = totDegree;

% ------------------------------------------------------------------------------
figure('color','w');
myColors = BF_getcmap('set1',5,1);

subplot(5,1,1); box('on'); hold on
% Degree distribution:
krAll = min(nodeData):max(nodeData);
N = arrayfun(@(x)sum(nodeData==x),krAll);
bar(krAll,N,'EdgeColor','k','FaceColor',myColors{4})
xlim([krAll(1)-0.5,krAll(end)+0.5]);
ylabel('Frequency')
xlabel('Degree')
title(sprintf('p threshold = %.2g',pThreshold),'interpreter','none')

% Rich curve:
subplot(5,1,2:5); box('on'); hold on
xlim([krAll(1)-0.5,krAll(end)+0.5]);
plot(totDegree,inDegree./totDegree,'xk')
xlabel('In + out degree')
ylabel('Proportion of in-links')

fprintf(1,'Average <k_in/k_tot> across nodes in the network: %f\n',mean(inDegree./totDegree));

top10Deg = quantile(totDegree,0.9);
isRich = (totDegree > top10Deg);
fprintf(1,'Average <k_in/k_tot> across top-10%% rich nodes: %f\n',mean(inDegree(isRich)./totDegree(isRich)));

set(gcf,'Position',[1000,756,539,582])