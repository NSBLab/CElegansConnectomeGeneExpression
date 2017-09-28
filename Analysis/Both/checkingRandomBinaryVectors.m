% function meanScore = checkingRandomBinaryVectors(n2)
N = 948;
n1 = 3;
n2 = 10;
nPerms = 100;
p = zeros(1,max([n1,n2]));

for i = 1:max([n1,n2])+1
    p(i) = probMMatches(N,n1,n2,i-1,true);
end
cump = cumsum([0,p]);

numMatches = zeros(nPerms,1);
scores = zeros(nPerms,1);
for i = 1:nPerms
    x1 = zeros(N,1);
    x1(randperm(N*1,n1)) = 1;
    x2 = zeros(N,1);
    x2(randperm(N*1,n2)) = 1;
    numMatches(i) = sum(x1==1 & x2==1);
    scores(i) = cump(numMatches(i)+1);
end
pEmpirical = arrayfun(@(x)mean(numMatches==x),0:max([n1,n2]));

figure('color','w'); ax = gca;
bar([p;pEmpirical]');
ax.XTickLabel = 0:max([n1,n2]);

cump1 = cumsum(p);
cump2 = cumsum([0,p(1:end-1)]);
cump = mean([cump1;cump2]);
meanScoreEmpirical = sum(cump.*pEmpirical)

% % Derive the theoretical mean score -- fewer than m matches:
% cump = cumsum([0,p]);
% cump = cump(1:end-1);
% meanScore1 = sum(cump.*p);
%
% % What about m matches or fewer:
% cump = cumsum(p);
% meanScore2 = sum(cump.*p);
%
% meanScore = mean([meanScore1,meanScore2]);

%% Mean score hybrid:
% cump1 = cumsum(p);
% cump2 = cumsum([0,p(1:end-1)]);
% cump = mean([cump1;cump2]);
% meanScore = sum(cump.*p);


% fprintf(1,'theoretical: %.3f, empirical: %.3f\n',meanScore,mean(scores));
% end


%===============================================================================
numChecks = 10;
N = 948;
n1 = 6;
n2 = 10;
for i = 1:numChecks
    x1 = zeros(N,1);
    x1(randperm(N*1,n1)) = 1;
    x2 = zeros(N,1);
    x2(randperm(N*1,n2)) = 1;

    X = [x1';x2'];

    coX = CoExpressionPositiveMatch(X,true,true);

end
