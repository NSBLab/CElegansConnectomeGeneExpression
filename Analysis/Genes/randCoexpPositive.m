
numRuns = 10; 
genes = G.GeneExpData.Direct; 
randgenes = zeros(size(genes,1), size(genes,2), numRuns);
for j=1:numRuns
    for i=1:size(genes,1)
        gene = genes(i,:);
        randgenes(i,:,j) = gene(randperm(length(gene)));
    end
end

coX = zeros(279,279,numRuns); 
for k=1:numRuns
    X = squeeze(randgenes(:,:,k)); 
    coX(:,:,k) = CoExpressionPositiveMatch(X,true); 
end

randCoexp = mean(coX,3); 
