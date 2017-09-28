function G = GeneName2EntrezID(G,whatData,annotationType)
% conversion made using https://david.ncifcrf.gov/conversion.jsp?VFROM=NA
% submitted list - 577 genes names from TEST.txt

fileName = 'CelegansEntrezID.txt';
fileID = fopen(fileName,'r'); % read in text file % Formatspec: if an error occurs for a different file, try regenerating the code from the Import Tool.
if fileID==-1
    error('Cannot read in %s (try regenerating the code from the Import Tool)',fileName);
end
fprintf(1,'Loaded entrez ID matching list from %s\n',fileName);

formatSpec = '%s%s%[^\n\r]';
delimiter = '\t';
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'ReturnOnError',false);
fclose(fileID);

dataArray{1} = dataArray{1}(2:end);
dataArray{2} = dataArray{2}(2:end);
dataArray{3} = dataArray{3}(2:end);

IDtable = [dataArray{1}, dataArray{2}];
% choose genes for C elegans, ignore others
k=1;
for i=1:length(dataArray{1})
    if strncmp('Caenorhabditis',dataArray{3}(i),14)
        ind(k)=i;
        k=k+1;
    end
end

IDtable = IDtable(ind,:);
geneNames = [G.geneAcronyms.Direct];
% for each gene in Gene structure input entrexID
for j=1:length(IDtable)
    Name1 = geneNames{j};
    for l=1:length(IDtable)
        if strcmpi(Name1,IDtable{l,1})
            match(j)=l;
        end
    end

end
% add entrez ID to Gene structure
for j=1:length(IDtable)
    IDs(j)=str2double(IDtable{j,2});
end

for p = 1:length(IDtable)
    G.GeneStruct(p).entrezID = IDs(match(p));
end
fileNameSave = sprintf('CelegansGeneDataWS%d%s.mat',whatData,annotationType);
fileNameSave = fullfile('Data',fileNameSave); % place in the Data directory
save(fileNameSave,'G','-append');
fprintf(1,'Entrez IDs for genes saved to %s\n',fileNameSave);

end
