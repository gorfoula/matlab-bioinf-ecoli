function [mutantfile]=SynthesizeMutantFile(MutantFile,AllSeqFile)

%% Read MUTANT
[text_] = textread(MutantFile,'%s',-1,'delimiter','\n'); % read all lines
header_mut=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(2:end),'[\t]'); % split coloumns
Table_MUT=CellTable2StrTable(Table_); % convert Cell table to string table
%% Read DataBase
[text_] = textread(AllSeqFile,'%s',-1,'delimiter','\n'); % read all lines
header_db=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(2:end),'[\t]'); % split coloumns
Table_DB=CellTable2StrTable(Table_); % convert Cell table to string table

% [Name, SP, Mature, FolwSeq, Efficiency, Ref, Comments] = textread(MutantFile,'%s %s %s %s %d %s %s',-1,'delimiter','\t');
% [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(AllSeqFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
% GnNames, Accessions, Comments, Descr, TotLen, MW, SPLen, SEQ

mutants=size(Table_MUT,1);
index=find(strcmp(Table_DB(:,1),'PPB_ECOLI'));
SEQUENCE=Table_DB(index,end);

AA_mut=Double2CellTable([1:1:mutants]);
[Names_mut]=StickEndStr(AA_mut','MUT','_');
Comments_mut=Table_MUT(:,1);
[Descr_mut]=MergeColumns([Table_MUT(:,end-1:end) Table_MUT(:,2)]);
TotLen_mut=zeros(mutants,1);
MW_mut=zeros(mutants,1);
SPLen_mut=zeros(mutants,1);
SEQ_mut=cell(mutants,1);

for i=1:1:mutants
%     i
%     Table_MUT{i,3}
    SPLen_mut(i)=length(Table_MUT{i,3});
    [s_pos e_pos]=regexp(SEQUENCE,Table_MUT{i,5});
    if(strcmp(Table_MUT{i,5},''))
        SEQ_mut{i}=[Table_MUT{i,3},Table_MUT{i,4},Table_MUT{i,5}];
    else
        SEQ_mut{i}=[Table_MUT{i,3},Table_MUT{i,4},Table_MUT{i,5},SEQUENCE{1}(e_pos{1}+1:end)];
    end
    
    MW_mut(i)=molweight(SEQ_mut{i});
    TotLen_mut(i)=length(SEQ_mut{i});
end
found_path=regexp(MutantFile,'[/]');
found_suffix=regexp(MutantFile,'[.]');
header='GnNames\tEfficiency\tDescr Names\tComments\tTotLen\tMW\tSPLen\tSEQ';
mutantfile=[MutantFile(1:found_path(end)),MutantFile(found_path(end)+1:found_suffix(end)-1),'_mn.txt']
FileWriteTable(mutantfile,[Names_mut,Table_MUT(:,end-2),Comments_mut, Descr_mut, Double2CellTable(TotLen_mut), Double2CellTable(MW_mut), Double2CellTable(SPLen_mut), SEQ_mut],header,'w');
FastaWrite(['Data/fASTA/',MutantFile(found_path(end)+1:found_suffix(end)-1),'.fasta'],Names_mut,SEQ_mut,70,'w');
% WriteMutantFiles(mutantfile,mutants,GnNames_mut, Descr_mut, TotLen_mut, MW_mut, SPLen_mut, SEQ_mut);

end

function [Comb_SEQ]=StickEndStr(SEQ,strL,strR)

[peptides]=size(SEQ,1);
Comb_SEQ=SEQ;
for i=1:peptides
    Comb_SEQ(i)={[strL,SEQ{i},strR]};
end

end