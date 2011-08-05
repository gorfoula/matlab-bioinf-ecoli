function []=CompareHydrophobicPatches(DatasetFiles,windows)
%% INPUT
%%  DatasetFiles:   File containing File Names of datasets to be compared
%%  windows:    lenght range of hydropho patches
%% INIT
close all;
[FileNames] = textread(DatasetFiles,'%s',-1,'delimiter','\t');
files=length(FileNames);
legend_names=cell(1,files);
len=120;
stpos=1:2:len;
numpatch=0:1:10;
patchlen=windows(1):windows(2);
pos=1:2:len;
Len_dist=zeros(length(patchlen),files);Len=[];
Start_Pos_dist=zeros(length(stpos),files);Start_Pos=[];
Num_Patches_dist=zeros(length(numpatch),files);Num_Patches=[];
Per_Pos_dist=zeros(length(pos),files);Per_Pos=[];
for i=1:1:files
    %% Locate Hydrophobic Patches
    [positions sequences patches AllPosPho]=HydrophobicPatches(FileNames{i},len,windows,'scale');
    proteins=length(patches);
    proteins_with_patches=sum(patches>0)/proteins; % percent of proteins that have at least one patch
    percentage(i)=proteins_with_patches*100;
    %% Legend Names
    [si ei ext matches]=regexp(FileNames{i},'[/_.]*[A-Za-z_0-9]+[.]');
    legend_names(i)={[matches{1}(2:end-1),' ',num2str(round(percentage(i))),'% have patches']};
    %% Figures
    
    [freq x]=hist(positions(:,1),patchlen);
    Len=[Len;[positions(:,1)  ones(length(positions(:,1)),1)*i]];
    Len_dist(:,i)=freq*proteins_with_patches*100./(length(positions(:,1)));  % Patches Length
    [freq x]=hist(positions(:,2),stpos);
    Start_Pos=[Start_Pos;[positions(:,2) ones(length(positions(:,2)),1)*i]];  
    Start_Pos_dist(:,i)=freq*proteins_with_patches*100./proteins;  % Starting position of patch
    [freq x]=hist(patches,numpatch);
    Num_Patches=[Num_Patches;[patches ones(length(patches),1)*i]];
    Num_Patches_dist(:,i)=(freq*100./proteins); % Number of patches
    [freq x]=hist(AllPosPho,pos);
    Per_Pos=[Per_Pos;[AllPosPho' ones(length(AllPosPho),1)*i]]
    Per_Pos_dist(:,i)=(freq*100./proteins);   % per position patch on mature
end

%% Figure titles and legends
[h1]=FigureLegends(patchlen,Len_dist,1,'Length (aas)','Percent','Width of Patches',legend_names,'p',{'-','';':','o'},Len(:,1),Len(:,2));
[h2]=FigureLegends(stpos,Start_Pos_dist,2,'Postion (aas)','Percent','Starting Position of Patches',legend_names,'pm',{'-','';':','o'},Start_Pos(:,1),Start_Pos(:,2));
[h3]=FigureLegends(numpatch,Num_Patches_dist,3,'Number of patches','Percent','Number of patches per Protein',legend_names,'b',{'-' 'o' },Num_Patches(:,1),Num_Patches(:,2));
[h4]=FigureLegends(pos,Per_Pos_dist,4,'Postion (aas)','Percent','Per Position residues',legend_names,'pm',{'-','';':','o'},Per_Pos(:,1),Per_Pos(:,2));

%% Save figures
saveas(h1,['Figures/Patches_',int2str(windows(1)),'-',int2str(windows(end)),'.bmp'],'bmp');
saveas(h2,['Figures/Paosition_',int2str(windows(1)),'-',int2str(windows(end)),'.bmp'],'bmp');
saveas(h3,['Figures/NoOfPatches_',int2str(windows(1)),'-',int2str(windows(end)),'.bmp'],'bmp');
saveas(h4,['Figures/PerPosPhoResidues',int2str(windows(1)),'-',int2str(windows(end)),'.bmp'],'bmp');

hold off;
end