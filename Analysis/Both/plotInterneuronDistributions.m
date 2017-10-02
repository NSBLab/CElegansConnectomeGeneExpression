function [dataCell, P,S, nC, nN] = plotInterneuronDistributions(C,G)

D = GiveMeDefault; 
kHub = D.kHub; 
maskLR = GiveMeLRMask(C); 
namesCommand = {'AVAL'; 'AVAR'; 'AVBL'; 'AVBR'; 'AVDL'; 'AVDR'; 'AVEL'; 'AVER'; 'PVCL'; 'PVCR'};
namesInterneurons = {'AVFL'; 'AVFR'; 'AVHL'; 'AVHR';'AVJL'; 'AVJR'; 'AVKR'}; 

Adj = GiveMeAdj(C, 'zeroBinary'); 
[~,~,deg] = degrees_dir(Adj); 

Names = C.NeuronNames; 

hubs = deg>kHub; 
%make masks for hubs
[~,commandnames] = intersect(Names, namesCommand); 
command = false(1,279); 
command(commandnames)=1; 
%make masks for non-hub interneurons
[~,interneuronnames] = intersect(Names, namesInterneurons); 
interneurons = false(1,279); 
interneurons(interneuronnames)=1; 

% exclude LR by making them 0.
coexp = GiveMeCoexpression(G,[],false).*~maskLR;


%figure; imagesc(COEXP); yticks([1:9,13]); yticklabels(Names(command)); 
%xticks([1:9,13]); xticklabels(Names(command)); 

%% coexpression for all hubs
allH = nonzeros(triu(coexp(hubs==1, hubs==1),1)); 

%% coexpression for matched non-hub interneurons
nonhubInterneurons = nonzeros(triu(coexp(interneurons==1, interneurons==1),1)); 

%% coexpression for hub command interneurons
hubANDcommand = hubs & command; 
nC = sum(hubANDcommand); 
commandH = nonzeros(triu(coexp(hubANDcommand==1, hubANDcommand==1),1)); 

%% coexpression between other hub neurons
hubNOTcommand = hubs & ~command; 
nN = sum(hubNOTcommand); 
noncommandH = nonzeros(triu(coexp(hubNOTcommand==1, hubNOTcommand==1),1)); 

dataCell{1} = allH; dataCell{2} = nonhubInterneurons; dataCell{3} = commandH; dataCell{4} = noncommandH; 
%JitteredParallelScatter(dataCell); 
[P.commandNONcommand, ~, S.commandNONcommand] = ranksum(commandH,noncommandH); 
[P.noncommandNONhub, ~, S.noncommandNONhub] = ranksum(noncommandH,nonhubInterneurons); 
end


