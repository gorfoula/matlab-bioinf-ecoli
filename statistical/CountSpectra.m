function []=CountSpectra(file,filedb,index)

%% This function given the output of scafold it calculates
%%      -> Unique Peptides
%%      -> Unique Spectra
%%      -> Total Number of spectra
%%
%% INPUT
%%      file:   path of output file of scafold, in tab delimited version
%%      filedb: database file
%%      index:  [H X Y Z1 Z2 A B C D]
%%              H:  how many lines are occupied by header
%%              X:  column number of Protein Accession (1)
%%              Y:  column number of Peptide sequence (2)
%%              Z1: Protein MW (3)
%%              Z2: column number of Peptide MW (4)
%%              A:  column number 1H counts (5)
%%              B:  column number 2H counts (6)
%%              C:  column number 3H counts (7)
%%              D:  column number 4H counts (8)
%%
%%  EXAMPLE: CountSpectra('Peptide Report for PRO008_1DG010_L1_BL21_allmods.txt','TableS4_BL21_DE3_ProteomeV8.txt');


if(exist('index','var')==0)
   index=[1 6 15 8 28 23 24 25 26]; % header | Accession | Peptide Sequence | Protein Mass 
end

[INTABLE] = ReadTable(file,'\n');
[namefile dir]=IsolateFileName({file});
outfile=[dir{1},namefile{1},'_NSAF.txt'];
INTABLE_noheader=INTABLE(index(1)+1:end,:);

[DBTABLE] = ReadTable(filedb,'\n');
DBTABLE=DBTABLE(:,[1 2 3 23 26]); %% entry name k12 | entry name BL21 | Description | Sequence

[TABLE Address]=Douplicates(file,index(1:2),INTABLE);

OUTTABLE=cell(max(Address),9);

for i=unique(Address)'
    cur_INTABLE=INTABLE_noheader(Address==i,index(2:end));
    
    OUTTABLE{i,1}=[cur_INTABLE{1,1},'_ECOBD'];  %% Protein Accession
    OUTTABLE{i,5}=cur_INTABLE{1,4};    %%  Protein Mass
    
    [found]=strcmp(DBTABLE(:,2),OUTTABLE{i,1});
    
    if(sum(found)>0)
        FOUNDTABLE=DBTABLE(found,:);
    else
        FOUNDTABLE=cell(1,5);
        FOUNDTABLE{:,2}=OUTTABLE{i,1};
    end
    
    OUTTABLE(i,1:4)=FOUNDTABLE(:,1:4);%%  Fished fields
    OUTTABLE{i,6}=length(FOUNDTABLE{:,5});  %%% Protein length
    
    SEQ=cur_INTABLE(:,2);
    [TABLE Address_seq]=Douplicates(file,[0 1],SEQ);
    OUTTABLE{i,7}=max(Address_seq);  %% number of unique peptides

    SEQ_MW_merged=MergeColumns(cur_INTABLE(:,[2,4]),'_');
    [TABLE Address_seqmwmerged]=Douplicates(file,[0 1],SEQ_MW_merged);
    uniq_modified_spectra=max(Address_seqmwmerged);  %% number of unique peptides
    
    Charges_count=CellTable2Double(cur_INTABLE(:,5:8));
    Charges_count_gt_zero=Charges_count>0;
    uniq_charges=sum(Charges_count_gt_zero,2);
    plus5charge=(uniq_charges==0);
    uniq_charges(uniq_charges==0)=1;
    
    uniq_spectra=0;
    for j=unique(Address_seqmwmerged)'
%         SEQ_MW_merged(Address_seqmwmerged==j)
%         Charges_count(Address_seqmwmerged==j,:)
%         uniq_charges(Address_seqmwmerged==j,:)

        uniq_spectra=uniq_spectra+max(uniq_charges(Address_seqmwmerged==j));
    end
    
    OUTTABLE{i,8}=uniq_spectra;  %% Number of unique spectra
    OUTTABLE{i,9}=sum(sum(Charges_count))+sum(plus5charge); %% Number of total spectra
    
end

FileWriteTable(outfile,[[{'Entry Name K12'} {'Entry Name BL21(DE3)'} {'Category'} {'Description'} {'Protein MW'} {'Protein Len'} {'Unique Peptides'} {'Unique Spectra'} {'Total #spectra'}];OUTTABLE],[],'w');

end