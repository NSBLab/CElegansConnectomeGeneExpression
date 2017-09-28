function connectionWeightVS(C,G, measureTOcompare)
Adj = GiveMeAdj(C,'zeroWeighted'); 
for i=1:279
    for j=1:279
        if Adj(i,j)~=0
            Adj(j,i) = Adj(i,j);
        end
    end
end

M = triu(logical(Adj),1);
if strcmp(measureTOcompare,'coexpression')
    measure = GiveMeCoexpression(G).*M;
elseif strcmp(measureTOcompare,'distance')
    measure = GiveMeDist(C).*M;
end



Adj = Adj.*M;
data1 =[measure(:), Adj(:)];
data1( ~any(data1,2), : ) = [];

figure;
if strcmp(measureTOcompare,'coexpression')
    scatter(data1(:,2), data1(:,1)); ylabel('Coexpression');xlabel('Connection weight');
    BF_PlotQuantiles(data1(:,2), data1(:,1),100,0,true); ylabel('Coexpression'); xlabel('Connection weight');
elseif strcmp(measureTOcompare,'distance')
   scatter(data1(:,2), data1(:,1)); ylabel('Distance');xlabel('Connection weight');
   BF_PlotQuantiles(data1(:,2), data1(:,1),100,0,true);ylabel('Distance');xlabel('Connection weight');
end

end
