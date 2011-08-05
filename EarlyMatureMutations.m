function []=EarlyMatureMutations(IFile,IFile_cyto,SelFeatFILE,WeigthFile,sel,area,mut,fig)

close all;

%%% Read Secreted
[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNames_cyto, Descr_cyto, TotLen_cyto, MW_cyto, SPLen_cyto, SEQ_cyto] = textread(IFile_cyto,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
%%% Read Selected Features File
[aa, FeatureID] = textread(SelFeatFILE,'%d %d',-1,'delimiter','\t');
FeatureID=FeatureID-1;  %%%% we don't want first column which is for category
NUMofFEATURES=aa(end);
%%%%% Read weight file
[w] = textread(WeigthFile,'%f',-1,'delimiter','\t');  %%% read weights
VALUES=20;  %%% division factor / number of attribute IDs
len=100;
%%% SECRETED Only Peptides with sufficient amino acid length
MOLECULARW=MW(TotLen>=SPLen+len);
DESCRIPTION=Descr(TotLen>=SPLen+len);
Names=GnNames(TotLen>=SPLen+len);   %Gene Names
SELECTED=SEQ(TotLen>=SPLen+len);    %Select Petides with Seq long enough
SIGNALS=length(SELECTED);
TOTALLEN=TotLen(TotLen>=SPLen+len);
SPLen=SPLen(TotLen>=SPLen+len);     %Select corresponding Signal Peptide Lengths
                                    % title of each column
%%%%  Open files for good and bad mutations
found=regexp(SelFeatFILE,'[/]');
switch (mut)
    case 0
        OutFile_good=['Data/GoodMutationsMat_CHOP_',SelFeatFILE(found(end)+1:end-4),'.txt'];
        OutFile_bad=['Data/BadMutationsMat_CHOP_',SelFeatFILE(found(end)+1:end-4),'.txt'];
    case 1
        OutFile_good=['Data/GoodMutationsMat_Shift_',SelFeatFILE(found(end)+1:end-4),'.txt'];
        OutFile_bad=['Data/BadMutationsMat_Shift_',SelFeatFILE(found(end)+1:end-4),'.txt'];
    otherwise
        display('=>Wrong type of mutation at <EarlyMatureMutations.m>');
        return;
end
s=fopen(OutFile_good,'w');
fprintf(s,'Name \t Mutant \t Shift or Chop \t Score \t SP \t Mature \n');
s=fopen(OutFile_bad,'w');
fprintf(s,'Name \t Mutant \t Shift or Chop \t Score \t SP \t Mature \n');
fclose all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step=5;
for i=1:1:SIGNALS
    display(Names{i});
    if(sel==20119)
       switch (mut)
           case 0
               [Combo]=ChopEarlyMature('',SELECTED{i}(SPLen(i)+1:end),step,len);x=(0:step:length(Combo)*step-step);  %%% chop only mature model
           case 1
               [Combo]=ShiftEarlyMature('',SELECTED{i}(SPLen(i)+1:len+SPLen(i)),area);x=1:1:length(Combo);   %%% Shift only mature model
           otherwise
            display('=>Wrong type of mutation at <EarlyMatureMutations.m>');
            return;
       end
       cur_SP='';
       cur_MATURE=SELECTED{i}(SPLen(i)+1:len+SPLen(i));
    else
        switch (mut)
           case 0
               [Combo]=ChopEarlyMature(SELECTED{i}(2:SPLen(i)),SELECTED{i}(SPLen(i)+1:end),step,len);x=0:step:length(Combo)*step-step;
           case 1
               [Combo]=ShiftEarlyMature(SELECTED{i}(2:SPLen(i)),SELECTED{i}(SPLen(i)+1:len+1),area);x=1:1:length(Combo);
           otherwise
            display('=>Wrong type of mutation at <EarlyMatureMutations.m>');
            return;
        end
        cur_SP=SELECTED{i}(1:SPLen(i));
        cur_MATURE=SELECTED{i}(SPLen(i)+1:len);
    end
    
    [allscores]=ScoreCalculation(Combo,w,FeatureID,VALUES,len,sel);   %%%% Score calculation
    %%%%% Print Mutation Score stem plot
    if(fig==1)   
        figure(1);
        stem(x,allscores);
        title(Names{i});
        pause;
    end
    %%%%% File write of worst and best mutations
    thres=(allscores>=3);
    FileWrite(OutFile_good,Names{i},Combo(thres),x(thres),allscores(thres),cur_SP,cur_MATURE);
    thres=(allscores<=0.5);
    FileWrite(OutFile_bad,Names{i},Combo(thres),x(thres),allscores(thres),cur_SP,cur_MATURE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end


function []=FileWrite(IFile,Name,Sequence,shift_chop,allscores,SP,MATURE)
s=fopen(IFile,'a');
combos=length(Sequence);
for j=1:1:combos   %for all best signal peptides
    fprintf(s,'%s \t %s \t %d \t %d \t %s \t %s \n',Name,Sequence{j},shift_chop(j),allscores(j),SP,MATURE);  %% file save best SP in an index        
end
fclose all;
end
