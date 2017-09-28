function coX = CoExpressionPositiveMatch(X,doExact,loadFile,talkToMe)
%-------------------------------------------------------------------------------
% Idea is to use positive matches for a pair of binary strings
% Binary matrix is X, repeats across all pairs of rows
% Output is coexpression matrix, coX
%-------------------------------------------------------------------------------
tic
if nargin < 2
    doExact = true;
    % Kind of necessary to get meaningful results
end
if nargin < 3
    loadFile = true;
end
if nargin < 4
    talkToMe = true;
end
%===============================================================================

numRows = size(X,1);

% Try to load from file
if loadFile
    L = load('preComputePositiveMatch.mat','probMatches','numOnes','numMatches','numGenes','doExact');
    fprintf('Pre-computed results loaded\n');
    if size(X,2)~=L.numGenes
        warning('%u doesn''t match the length of vector stored in the .mat file (%u)',...
                                    size(X,2),L.numGenes)
        clear('L')
    elseif doExact~=L.doExact
        warning('loaded file doesn''t match doExact -- cannot load data')
        clear('L')
    end
end

% Compare every pair of rows into pairwise similarity matrix, coX
coX = zeros(numRows);
for i = 1:numRows-1
    for j = i:numRows
        coX(i,j) = permMatchScore(X(i,:),X(j,:));
    end
end
% Symmetrize:
coxD = diag(coX);
coX = coX + coX';
% Add diagonal back:
coX(logical(eye(size(coX)))) = coxD;

%-------------------------------------------------------------------------------
    function pScore = permMatchScore(x1,x2)
        % Compute number of matches:
        numMatches = sum(x1 & x2);
        numOnes = [sum(x1),sum(x2)];
        n1 = min(numOnes);
        n2 = max(numOnes);
        N = length(x1);

        % Given binary vectors x1 and x2, computes probability of matches
        ps = zeros(numMatches+1,1);
        for ii = 1:numMatches+1
            if exist('L','var') & ~isnan(L.probMatches(L.numOnes==n1,L.numOnes==n2,L.numMatches==ii-1))
                if talkToMe
                    fprintf(1,'Loaded result from file for (%u,%u/%u) %u matches\n',n1,n2,N,ii-1);
                end
                ps(ii) = L.probMatches(L.numOnes==n1,L.numOnes==n2,L.numMatches==ii-1);
            else
                if talkToMe
                    fprintf(1,'Computing for (%u,%u/%u) %u matches\n',n1,n2,N,ii-1);
                end
                ps(ii) = probMMatches(N,n1,n2,ii-1,doExact);
            end
        end
        % pScore is the probability of getting fewer matches than that observed, by chance
        pScore = sum(ps(1:numMatches))+0.5*ps(end);
        % p = factorial(n1)*factorial(n2)*factorial(N-n1)*factorial(N-n2)/...
        %         (factorial(N)*factorial(n2-numMatches)*factorial(numMatches)*...
        %             factorial(N-n1-n2+numMatches));
    end
toc
end
