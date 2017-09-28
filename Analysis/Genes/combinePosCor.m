% red in positive coexpression measure matrices and make average
cd ('/projects/kg98/aurina/correlation')
coexp = zeros(100,150,150); 
for i=1:100
    load(sprintf('%d.mat', i)); 
    coexp(i,:,:) = BinCorrelation;
end

mcoexp = squeeze(mean(coexp,1)); 
figure; imagesc(mcoexp); 