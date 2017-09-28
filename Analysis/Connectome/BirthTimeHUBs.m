% compare birth times for feed in nodes (non-hubs sending inf to hubs), feed-out nodes (non-hubs receiving inf from hubs), hubs and non-hubs at each degree threshold.
% birth time
Adj = GiveMeAdj(C,'zeroBinary'); 
[kin,kout,deg] = degrees_dir(Adj);
birthTime = C.BirthTime;

BHub = zeros(max(deg),1);
BNONhub = zeros(max(deg),1);
Bfeedin = zeros(max(deg),1);
Bfeedout = zeros(max(deg),1);

%% calculate average
for k=1:max(deg)
isHub = deg>=k;
% mean birth time for hubs and non-hubs at each threshold
BHub(k) = mean(birthTime(isHub));
BNONhub(k) = mean(birthTime(~isHub));
% mean birth time for feed-in neurons and feed-out neurons
%find all neurons that feed-in to hubs
feedin = logical(sum(Adj(:,isHub==1),2));
% calculate 'weighted average' of birth times
Bfeedin(k)=mean(nonzeros(birthTime.*feedin));
%find all neurons that feed-out to hubs
feedout = logical(sum(Adj(isHub==1,:),1));
% calculate 'weighted average' of birth times
Bfeedout(k)=mean(nonzeros(birthTime.*feedout'));
end

figure;
scatter(1:max(deg), BHub,  'MarkerEdgeColor' ,[1 .15 0.07]); hold on;
scatter(1:max(deg), BNONhub,'MarkerEdgeColor' , [0 .14 0.53]); hold on;
scatter(1:max(deg), Bfeedin, 'MarkerEdgeColor' , [1 .65 0]); hold on;
scatter(1:max(deg), Bfeedout, 'MarkerEdgeColor', [.47 .87 .47]);
xlabel('Hub threshold'); ylabel('Average birth time (min)');
legend('hubs','non-hubs)', 'feed-in', 'feed-out');
% calculate and compare distributions for feed-in and feed-out neurons at
% each hub threshold
zVal = zeros(max(deg),1);
pVal = zeros(max(deg),1);

for k=1:max(deg)
isHub = deg>=k;
% birth time for feed-in neurons and feed-out neurons
% find all neurons that feed-in to hubs
feedin = logical(sum(Adj(:,isHub==1),2));
% calculate 'weighted average' of birth times
Bfeedin=nonzeros(birthTime.*feedin);
%find all neurons that feed-out to hubs
feedout = logical(sum(Adj(isHub==1,:),1));
% calculate 'weighted average' of birth times
Bfeedout=nonzeros(birthTime.*feedout');
BH = birthTime.*isHub';
BNH = birthTime.*~isHub';
dataCell{1} = nonzeros(Bfeedin);
dataCell{2} = nonzeros(Bfeedout);
dataCell{3} = nonzeros(BH);
dataCell{4} = nonzeros(BNH);
%JitteredParallelScatter(dataCell);
[p,h,stats] = ranksum(dataCell{1}, dataCell{2});
zVal(k) = stats.zval;
pVal(k)=p;
end
figure; scatter(1:max(deg), zVal); hold on; scatter(1:max(deg), pVal); hold on;
plot([1 max(deg)] , [1 1]*0.05, 'LineWidth', 1, 'Color', 'r');
xlabel('Hub threshold'); ylabel('Zvalue and Pvalue');
title('Comparing birth times between feed-in and feed-out neurons');
legend('Zval comparing feed-in vs feed-out birth times','Pval');





%end
