function RegenerateData()
%-------------------------------------------------------------------------------
% Idea is to read in all of the raw data files, process and organize them, then
% save them to the necessary output files
%-------------------------------------------------------------------------------

% First generate connectivity information -> save to file (for C structure)
fprintf(1,'\n\n\n~~FIRST CONNECTIVITY INFORMATION~~~\n\n\n');
Print_neuronconnect();

% Then use to match to gene expression information -> save to file (for G structure)
fprintf(1,'\n\n\n~~NOW TO PROCESS GENE EXPRESSION INFORMATION~~~\n\n\n');
SaveExpressionData();

end
