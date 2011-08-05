clear all;
close all;

Proteins_per_Cell=[14.22 660.47 2693.30 1220.90 381.36 22.71 93.79 11.36 783.08 4990.02 3414.65 20.45 639.52 790.41 31518.82 1350.98 128.00 679.82 1353.59 22075.01 9.44 1842.38];

f1='JoinFiles/Peripheral/30S_sel.txt';
f2='JoinFiles/Peripheral/UniprotRibosomalUnits.sif';
[FOUNDTABLE]=FishLines(f1,f2,[0 1],[1 3],-1,1);


f1='JoinFiles/GOAnnotations/ORList_Merge_GOs.txt';
f2='JoinFiles/Peripheral/TechnicalR_all_v9_gorfo.txt';
[FOUNDTABLE]=FishLines(f1,f2,[1 1],[1 2],-1,1);
JoinMSMSResults('JoinFiles/Peripheral/Merge_CellularComp_Biological.txt',1);


f1='JoinFiles/BL2119/TableS7_BL21_DE3_ProteomeV10.txt';
f2='JoinFiles/k12/TableS3_K12_ProteomeV18.txt';
[FOUNDTABLE]=FishLines(f2,f1,[1 1],[1 1],5,1);

f1='JoinFiles/k12/TableS3_K12_ProteomeV15.txt';
f2='JoinFiles/Peripheral/TechnicalR_all_v9_gorfo.txt';
f3='JoinFiles/Peripheral/Merge_CellularComp_Biological.txt';
[FOUNDTABLE]=FishLines(f1,f2,[1 1],[1 1],2,1);
[FOUNDTABLE]=FishLines(f3,f2,[1 1],[1 3],-1,1);


f1='JoinFiles/Peripheral/mistake.txt';
f2='JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt';
[FOUNDTABLE]=FishLines(f1,f2,[0 1],[1 2],1,1);

f2='JoinFiles/Peripheral/TechnicalR_all_v7.txt';
f1='JoinFiles/Peripheral/ProteinsPerMolecules.txt';
[FOUNDTABLE]=FishLines(f2,f1,[0 1],[1 1;2 2],-1,1);

file1='JoinFiles/k12/TableS3_K12_ProteomeV15.txt';
file2='JoinFiles/Peripheral/TechnicalR_all_v7.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 1],[1 1],[2],1);


file1='JoinFiles/k12/TableS3_K12_ProteomeV15.txt';
file2='JoinFiles/Peripheral/TechnicalR_all_v5_eco.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 1],[1 1],[2],1);

file1='JoinFiles/Peripheral/TechnicalR_all_v5_eco_UniprotAcc.txt';
file2='JoinFiles/Peripheral/Combination_EcoCyc_EcoWiki_Uniprot.tsv';
[FOUNDTABLE]=FishLines(file1,file2,[1 1],[1 2],[2 3],1);

file1='JoinFiles/Peripheral/GO_IDs_translation.txt';
file2='JoinFiles/Peripheral/TechnicalR_all_v5_eco_UniprotAcc_GO.txt';
[FOUNDTABLE]=FishLines(file1,file2,[0 1],[1 7],-1,1);

file1='JoinFiles/Peripheral/TechnicalR_all_v5_eco.txt';
% file2='JoinFiles/Peripheral/TechnicalR_translation.txt';
file2='JoinFiles/Peripheral/TechnicalR_all_v5_eco_BL21_GO_oneline.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 0],[2 1],-1,1);


hist(Proteins_per_Cell);
xlabel('Proteins Per molecule');
ylabel('Proteins');

file1='JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt';
file2='JoinFiles/Transcriptomics/Wei2001_LB_total.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 1],[2 1;3 1;4 1],-1,1);


file1='JoinFiles/k12/TableS3_K12_ProteomeV15.txt';
file2='JoinFiles/Peripheral/72InnerMem_NoTMsPredicted_TassosChecked.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 1],[1 1],-1,1);

file1='JoinFiles/k12/TableS3_K12_ProteomeV14.txt';
file2='JoinFiles/BernselVSCurated/IM_theoretical_Bernsel.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 0],[2 1],-1,1);

file1='JoinFiles/Peripheral/TechnicalR_all.txt';
file2='JoinFiles/Peripheral/InnerMem_NoTMs.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 0],[1 2],-1,1);

file1='JoinFiles/Peripheral/InnerMem_NoTMs.txt';
file2='JoinFiles/k12/TableS3_K12_ProteomeV14.txt';
[FOUNDTABLE]=FishLines(file2,file1,[1 0],[1 2],[1:29],1);

file1='JoinFiles/Peripheral/DetectedPeripheral.txt';
file2='JoinFiles/Peripheral/InnerMem_NoTMs.txt';
[FOUNDTABLE]=FishLines(file1,file2,[1 0],[1 2],[1:23],1);

file='JoinFiles/Peripheral/RibosomalIdentified.txt';
[FOUNDTABLE]=FishLines('JoinFiles/Peripheral/ProteinsPerMolecules.txt',file,[0 0],[1 1],[3 4 5 6 7 8],1);

file='JoinFiles/Peripheral/List_potential_F1.txt';
[FOUNDTABLE]=FishLines('JoinFiles/Peripheral/ProteinsPerMolecules.txt',file,[0 1],[1 1;2 2],[3 4 5 6 7 8],1);

file='JoinFiles/Peripheral/List_potential_F1.txt';
[FOUNDTABLE]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 1],[11],1);
file='JoinFiles/Peripheral/ProteinsPerMolecules.txt';
[FOUNDTABLE]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 0],[1 1;2 2],[11],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileY='JoinFiles/WTvsSecy/Proteins_1DG010_SecY.txt';
fileWT='JoinFiles/WTvsSecy/Proteins_1DG011_wt.txt';
[FOUNDTABLE]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',fileY,[1 1],[1 1],[2],1);
[FOUNDTABLE]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',fileWT,[1 1],[1 1],[2],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file='JoinFiles/k12/Spectral counting.txt';
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 1],[2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Merge transciptomic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  proteomics
TranscrDistSript('JoinFiles/Transcriptomics/Merge_LB.txt','JoinFiles/Transcriptomics/Merge_Protein_LB.txt','JoinFiles/Transcriptomics/Merge_mRNA_LB.txt',4)

file='JoinFiles/k12/Intact_83333_NoTags.txt';
[temp]=FishLines(file,'JoinFiles/k12/TableS3_K12_ProteomeV13_BasicProteome_3900.txt',[1 0],[1 2;3 2],[],1);

file='JoinFiles/k12/Peripheral.txt';
Proteome=ReadTable('JoinFiles/k12/TableS3_K12_ProteomeV13.txt');

[COMMON]=FishMergeInteractions([{'JoinFiles/Peripheral/Ecoli20101010_dip_interact_UPDATED_.txt'} {[1 1 2]}],[{'JoinFiles/k12/Cyto.txt'} {[1 2]}],[{'JoinFiles/k12/Peripheral.txt'} {[1 2]}]);
[COMMON_Entry]=MatchAccessionsToEntryName([{Proteome} {[1 2 1 5]}],Interactions);
FileWriteTable('JoinFiles/Peripheral/F1_A_PPI.txt',COMMON_Entry,[],'w');

[Interactions Names]=FishInteractions([{'JoinFiles/Peripheral/Ecoli20101010_dip_interact_UPDATED_.txt'} {[1 1 2]}],[{file} {[0 2]}]);
[COMMON_Entry]=MatchAccessionsToEntryName([{Proteome} {[1 2 1 5]}],Interactions);
FileWriteTable('JoinFiles/Peripheral/F1_PPI.txt',COMMON_Entry,[],'w');
Table=ReadTable(file);
[Table_Entry]=MatchAccessionsToEntryName([{Proteome} {[1 2 1 5]}],Table(:,2));
FileWriteTable('JoinFiles/Peripheral/F1_PPI_SList.txt',Table_Entry,[],'w');
%%%%%%%%%% Interactions  with A or F1 %%%%%%%%%%%%%%%%%%%
file='JoinFiles/Peripheral/Traped_Peripheral.txt';
filedb='JoinFiles/Peripheral/Cyto_Peripheral_Interactions.txt';
[temp]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 0],[1 1],[2]);
file='JoinFiles/Peripheral/Traped_Peripheral_Accession.txt';


file1='JoinFiles/Michalis/ColB.txt';
file2='JoinFiles/Michalis/ColA.txt';
[FOUNDTABLE11]=FishLines(file1,file2,[1 1],[1 1],[]);

[FOUNDTABLE11]=FishLines(file,filedb,[1 1],[1 1;1 2],[1]);
FileWriteTable('JoinFiles/Peripheral/Traped_Peripheral_Accession_Interactions.txt',[FOUNDTABLE11],[],'w');

file='JoinFiles/Peripheral/Traped_Peripheral_Accession_Interactions.txt';
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[2 1],[1 5]);
[FOUNDTABLE12]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[2 2],[1 5]);
%%%%%%%%%%% Fish Interactions of Interest   %%%%%%%%%%%%%
file='JoinFiles/k12/Cyto_OR_Peripheral.txt'; %%%% A and F1
[FOUNDTABLE11]=FishLines(file,'JoinFiles/Peripheral/Ecoli20101010_dip_interact_UPDATED_.txt',[0 1],[2 1;2 2],[1]);
FileWriteTable('JoinFiles/Peripheral/Cyto_Peripheral_Interactions.txt',[FOUNDTABLE11],[],'w');
file='JoinFiles/Peripheral/Potential_Peripheral.txt';  %%%% Potential F1
[temp]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 1],[2]);
file='JoinFiles/Peripheral/Potential_Peripheral_Accession.txt';
[FOUNDTABLE11]=FishLines(file,'JoinFiles/Peripheral/Ecoli20101010_dip_interact_UPDATED_.txt',[1 1],[1 1;1 2],[1]);
FileWriteTable('JoinFiles/Peripheral/Peripheral_Interactions.txt',[FOUNDTABLE11],[],'w');
%%%% Output: F1_Both_Potential_and_Curated.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  Fish Entry Names  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Not found Accessions are discarded  %%%%%%%%
file='JoinFiles/Peripheral/F1_Both_Potential_and_Curated.txt';
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[2 1],[1 5]);
[FOUNDTABLE12]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[2 2],[1 5]);
%%%% Output: F1_Both_Potential_and_Curated_ENTRYNAME.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Merge Entry Name with Curated Sub Loci %%%%%%%%%
file='JoinFiles/Peripheral/F1_Both_Potential_and_Curated_ENTRYNAME.txt';
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 1],[5]);
[FOUNDTABLE12]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 2],[5]);
New_Label=MergeColumns([FOUNDTABLE11(:,2) FOUNDTABLE11(:,1)],'_');
New_Label2=MergeColumns([FOUNDTABLE12(:,3) FOUNDTABLE12(:,1)],'_');
FileWriteTable('JoinFiles/Peripheral/F1_Both_Potential_and_Curated_ENTRYNAME_SUBLOCI.txt',[New_Label New_Label2 FOUNDTABLE11(:,2:end)],[],'w');
%%%% Output: F1_Both_Potential_and_Curated_ENTRYNAME_SUBLOCI.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Peripheral Entry Name Loci %%%%%%%%%%%%%
file='JoinFiles/Peripheral/Potential_Peripheral_Accession.txt';
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',file,[1 1],[1 2],[5]);
New_Label=MergeColumns([FOUNDTABLE11(:,3) FOUNDTABLE11(:,1)],'_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/AllPeri_AND_Pot_Peri_Interactions.txt',[1 1],[1 1],[5]);
[FOUNDTABLE12]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/AllPeri_AND_Pot_Peri_Interactions.txt',[1 1],[1 2],[5]);

New_Label=MergeColumns([FOUNDTABLE11(:,2) FOUNDTABLE11(:,1)],'_');
New_Label2=MergeColumns([FOUNDTABLE12(:,3) FOUNDTABLE12(:,1)],'_');
FileWriteTable('JoinFiles/Peripheral/AllPeri_AND_Pot_Peri_Interactions_SUBLOCI.txt',[New_Label New_Label2 FOUNDTABLE11(:,2:end)],[],'w');

[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Potential_Peripheral.txt',[1 1],[1 1],[5]);

[FOUNDTABLE12]=FishLines('JoinFiles/Peripheral/Potential_Peripheral.txt','JoinFiles/Peripheral/TableS3_K12_ProteomeV13(for)Cyto_Peripheral_Interactions.txt',[1 1],[1 1;1 2],[1 2 3 4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[FOUNDTABLE11]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Cyto_Peripheral_Interactions.txt',[1 1],[2 2],[1]);
[FOUNDTABLE12]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Cyto_Peripheral_Interactions.txt',[1 1],[2 3],[1]);

file1='JoinFiles/Peripheral/AllAccessions_from_Interactions_RetrieveUniprotAccession.txt';
file2='JoinFiles/Peripheral/Ecoli20101010_dip_interact.txt';
[FOUNDTABLE1]=FishLines(file1,file2,[1 1],[1 1]);
[FOUNDTABLE2]=FishLines(file1,file2,[1 1],[1 2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file='JoinFiles/Peripheral/Ecoli20101010_dip_interact.txt';

file2='JoinFiles/Peripheral/AllAccessions_fro_Interactions.txt';

[namefile dir]=IsolateFileName({file2});
outfile=[dir{1},namefile{1},'_Unique.txt'];
outfile2=[dir{1},namefile{1},'_Interactions.txt'];

[CATCH] = ReadTable(file2,'\n');
[TABLE Address]=Douplicates(file2,[1 1],CATCH);
FileWriteTable(outfile2,TABLE,'Uniprot Accession\t#Interactions','a');
CorespLink=MergeColumns(TABLE(:,1),'http://www.uniprot.org/uniprot/?query=replaces:','h');
CorespLink=MergeColumns(CorespLink,'+AND+ECOLI','t');
FileWriteTable(outfile,CorespLink,[],'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%
FastaSubSelection('JoinFiles/Peripheral/InnerMem_NoTMs.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[2 0]);
JoinMSMSResults('JoinFiles/k12/temp.txt',1);
%%%%%%%%%%%%%%%%%%%%%%%%
FishLines('JoinFiles/Peripheral/(DB)TableS3_K12_ProteomeV13(for)Potential_Peripheral(1)F.txt','JoinFiles/Peripheral/(DB)TableS3_K12_ProteomeV13(for)Ecoli20101010_dip_interact_updated(2)F.txt',[1 1;1 2],[1 1]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Ecoli20101010_dip_interact_updated.txt',[1 1],[1 1]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Potential_Peripheral.txt',[1 1],[1 1]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Peripheral/Ecoli20101010_dip_interact.txt',[2 1],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
SignalPeptideMatureFusions({'Gems/Curated/AllSecr.txt' 'Gems/Curated/AllCyto.txt'},1);
%%%%%%%%%%%%%%%%
[INTABLE] = ReadTable('JoinFiles/Files/PhoA_Mature_combos.txt','\n');
SCORES=CellTable2Double(INTABLE(:,2));
% SCORES=(SCORES./max(abs(SCORES)))+1;

indx=strcmp(INTABLE(:,1),'PPB_ECOLI');
phoA_score=SCORES(find(indx));

x=(-1:0.1:1)*max(abs(SCORES));

[freq]=hist(SCORES,x);
bar(x,(freq*100)./sum(freq));

x_phoa=ones(length(x),1)*phoA_score;
y_phoa=ones(length(x),1)*((freq*100)./sum(freq));

hold on; plot(x_phoa,y_phoa,'.-r')


xlabel('Secretion Score');
ylabel('Percent over all combinations');



%%%%%%%%%%%%%%%%%
JoinExperiments('JoinFiles/Files/Nikos_Fractions_1to24.txt',[1 1]); %%%% for Nikos Unique Peptides, fractions
JoinExperiments('JoinFiles/Files/Pro008_OG01_OG02_Fractions.txt',[1 1]); %%%% for Nikos Unique Peptides, fractions
%%%%%%%%%%%%%%%%%%
JoinMSMSResults('JoinFiles/GOAnnotations/Merge_GOs.txt',1);%% out file

TechnicalR_all_v9_gorfo.txt

JoinMSMSResults('JoinFiles/Files/Merge_Lanes.txt',1);%% out file
ExludeCyto('JoinFiles/Files/Merge_Lanes.txt',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SelectedFeatures('Gems/Curated/CytoVSsecr_FIN.txt',0,'Gems/Curated/TRAIN/CytoVSsecr.gemsin','Gems/Curated/TEST/CytoVSsecr.test',0,1);

SelectedFeatures('JoinFiles/Peripheral/Hiton_60_poly.txt',0,'JoinFiles/Peripheral/Peripheral_March2011_Cytoplasmic_March2011_gemsin.txt','JoinFiles/Peripheral/Peripheral_March2011_Cytoplasmic_March2011_gemsin.txt',0,1);
SelectedFeatures('JoinFiles/Peripheral/Hiton_100_poly.txt',0,'JoinFiles/Peripheral/Peripheral_March2011_Cytoplasmic_March2011_gemsin_100.txt','JoinFiles/Peripheral/Peripheral_March2011_Cytoplasmic_March2011_gemsin_100.txt',0,1);
%%%%%%%%%%%%%%%%%%%%%%
FastaSubSelection('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[1 1]);
FastaSubSelection('JoinFiles/k12/Peripheral.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[1 1]);
FastaSubSelection('JoinFiles/k12/InnerMembr.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[1 1]);
FastaSubSelection('JoinFiles/k12/Cyto.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HEAD SEQ]=FastaRead('JoinFiles/BL2119/uniprot-organism_469008.fasta');
FileWriteTable('JoinFiles/BL2119/FishSequence_fromfasta.txt',[[{'Entry Name K12'} {'Sequence'}];[HEAD' SEQ']],[],'w');
FishLines('JoinFiles/BL2119/FishSequence_fromfasta.txt','JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt',[1 2],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CountSpectra('JoinFiles/Files/Peptide Report for PRO008_1DG010_L1_BL21_allmods.txt',[1 6 15 8 28 24 25 26]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file='JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only(1)F(1)F.txt';
PrinteIDOnce(file,[1 1 2 3 4 11 12]);
file='JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides(1)F(1)F.txt';
PrinteIDOnce(file,[1 1 2 3 4 11 12]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Douplicates('JoinFiles/Total proc level1_ PRO008_1DG010_L1_BL21.txt',[1 1]);
FishLines('JoinFiles/Transcriptomics/TEMP1.txt','JoinFiles/Transcriptomics/TEMP2.txt',[1 1],[0 0]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13_BasicProteome_3900.txt','JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt',[1 1],[1 1]);

CompareExpTheoPeptides('JoinFiles/TheoreticalPeptides/Theoretical.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel(1)L.txt',[1 2 1 2 1 2],[1 1 1],1);
CompareExpTheoPeptides('JoinFiles/TheoreticalPeptides/Theoretical.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel(1)F.txt',[1 2 1 2 1 2],[1 1 1],1);

CompareExpTheoPeptides('JoinFiles/Transcriptomics/DetByUs_Cyto(1)L.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel(1)L.txt',[1 2 1 2 1 2],[1 1 1],1);
CompareExpTheoPeptides('JoinFiles/Transcriptomics/DetByUs_Cyto(1)L.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel(1)F.txt',[1 2 1 2 1 2],[1 1 1],1);

CompareExpTheoPeptides('JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Common_seq.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel(1)L.txt',[1 2 1 2 1 2],[1 1 1],1);
CompareExpTheoPeptides('JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Common_seq.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel(1)F.txt',[1 2 1 2 1 2],[1 1 1],1);

CompareExpTheoPeptides('JoinFiles/Transcriptomics/DetByUs_CEP(1)L.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel(1)L.txt',[1 2 1 2 1 2],[1 1 1],1);
CompareExpTheoPeptides('JoinFiles/Transcriptomics/DetByUs_CEP(1)L.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel(1)F.txt',[1 2 1 2 1 2],[1 1 1],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FishLines('JoinFiles/BL2119/uniprot-organism_469008.txt','JoinFiles/BL2119/temp.txt',[6 2],[1 1]);
FishLines('JoinFiles/Mobile Elements/911_HT.txt','JoinFiles/Mobile Elements/1006_HT.txt',[1 1],[1 1]);

FishLines('JoinFiles/CompareTrypticPeptides/BL21_Det_CEP.txt','JoinFiles/CompareTrypticPeptides/BL21_CEP.txt',[4 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Mobile Elements/EcoGene_IS_Prophage.txt',[2 1;3 1;4 6;1 5],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/NotDetByUs_Morethan10ProteinCopies.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/CEP_NonDetected_Common.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/BL21_CEP_Common.txt',[6 1],[1 1]);


FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Mobile Elements/EcoGene_IS_Prophage(2)F.txt',[2 1],[1 1]);

FishLines('JoinFiles/Mobile Elements/K12_Mixed_info_DB.txt','JoinFiles/Mobile Elements/700_Mobile_Elements_S1.txt',[2 2;3 2;4 3;12 5;13 6],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Mobile Elements/700_Mobile_Elements_S1.txt',[2 2;3 2;4 3],[1 1]);
FishLines('JoinFiles/Mobile Elements/EcoData022711-074325.txt','JoinFiles/Mobile Elements/700_Mobile_Elements_S1.txt',[2 2;6 5;7 6],[1 1]);
FishLines('JoinFiles/Mobile Elements/EcoData022711-074325.txt','JoinFiles/Mobile Elements/700_Mobile_Elements_S1(2)F.txt',[2 2;2 3;1 3;5 4;6 5;7 6],[1 1]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/Mobile Elements/700_Mobile_Elements_S1(2)F.txt',[2 1],[1 1]);
FishLines('JoinFiles/Mobile Elements/700_Mobile_Elements_S1(2)F.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 2],[1 1]);

FishLines('JoinFiles/Mobile Elements/K12_Mobile(2)F(for)TableS3_K12_ProteomeV13(1)F.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 2],[1 1]);
FishLines('JoinFiles/Mobile Elements/MobileElements.txt','JoinFiles/k12/TEMP.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/CEP.txt',[1 1],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/CYTO.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Transcriptomics/Non_DetByUs_CEP.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_Matched.txt',[2 1],[1 1]);

FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides.txt',[2 1],[1 1]);

FishLines('JoinFiles/ELOCATION/(DB)uniprot-organism_83333_SubLocFixed(for)CoBaltTable2(2)F.txt','JoinFiles/k12/TEMP.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/ELOCATION/CoBaltTable2.txt',[2 4;3 4],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[2 1],[1 1]);

FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/TheoreticalPeptides/Theoretical_CHECK_IDs_Theo(1)F.txt',[2 1],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/TheoreticalPeptides/Theoretical_CHECK_IDs_Theo.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/DetByUs_Cyto.txt',[1 1],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/DetByUs_CEP.txt',[1 1],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel.txt',[1 1;],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel.txt',[1 1;],[1 0]);

[list,Header,files,discr]=ExludeCyto('JoinFiles/OtherStudies/Notepad Lists/InFiles_All studies.txt',36);
[list,Header,files,discr]=ExludeCyto('JoinFiles/OtherStudies/Notepad Lists/Studies_wide_10.txt',13);

ExludeCyto('JoinFiles/OtherStudies/Notepad Lists/Wide_8(over_500)_LB.txt',13);


[CELLENV_wide_studies files_wide_studies CYTO_studies]=ExludeCyto('JoinFiles/OtherStudies/Notepad Lists/Merge_BL21_only.txt',2,1);

ExpressionAnalysis('JoinFiles/Transcriptomics/Merge_Protein_LB.txt','JoinFiles/Transcriptomics/Merge_mRNA_LB.txt',1);

[CELLENV_wide_studies files_wide_studies CYTO_studies]=ExludeCyto('JoinFiles/Transcriptomics/Merge_Wide_VS_Focused.txt',3,1);
%%%%%%%%%%%  COMPARE PROTEOMICS  %%%%%%%%%

JoinMSMSResults('JoinFiles/CompareTrypticPeptides/Merge_temp.txt',1);%% out file

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/CompareTrypticPeptides/ORList_Merge_temp.txt',[6 1],[1 1]);

CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/CompareTrypticPeptides/(DB)uniprot-organism_83333(for)ORList_Merge_temp(6)F.txt',[],6,[1 1],[1 1]);
CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/CompareTrypticPeptides/LH_Proteins.txt',[],6,[1 1],[1 0]);
CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/k12/CatgID.TXT','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[],5,[1 1],[1 1]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JoinMSMSResults('JoinFiles/Transcriptomics/Merge_Protein_LB_measurements.txt',1);%% out file
ExludeCyto('JoinFiles/Transcriptomics/Merge_Wide.txt',3,1);%% file ref | rerun
FishLines('JoinFiles/Transcriptomics/ORList_Merge_Protein_LB_measurements.txt','JoinFiles/Transcriptomics/(Merge_Wide)_COMPARISON_CELL_ENV_MARKED.txt',[1 1],[1 4]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/Transcriptomics/(DB)ORList_Merge_Protein_LB_measurements(for)(Merge_Wide)_COMPARISON_CELL_ENV_MARKED(1)F.txt',[1 12],[1 1]);

ExludeCyto('JoinFiles/OtherStudies/Notepad Lists/OthersFiles.txt',38,1);%% file ref | rerun
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/OtherStudies/Notepad Lists/(OthersFiles)_COMPARISON_CELL_ENV_MARKED.txt',[1 41],[1 4]);
% ExludeCyto('JoinFiles/Transcriptomics/Merge_Focused.txt',3,1);%% file ref | rerun
% FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/Transcriptomics/(Merge_Focused)_COMPARISON_CELL_ENV_MARKED.txt',[1 21],[1 4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[list,Header,files,discr]=JoinMSMSResults('JoinFiles/OtherStudies/Medium/LB_Merge.txt',1);
[list,Header,files,discr]=JoinMSMSResults('JoinFiles/OtherStudies/Notepad Lists/Wide_8(over_500)_LB.txt',1);
[list,Header,files,discr]=JoinMSMSResults('JoinFiles/OtherStudies/Notepad Lists/Focused_15(below_500)_LB.txt',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WideStudies='JoinFiles/OtherStudies/Notepad Lists/Studies_wide_10(over_500).txt';
[CELLENV_wide_studies files_wide_studies CYTO_studies]=ExludeCyto(WideStudies,12,1);
FishLines('JoinFiles/Transcriptomics/(Merge_Protein_LB)_COMPARISON_CELL_ENV_MARKED.txt','JoinFiles/OtherStudies/Notepad Lists/(Studies_wide_10(over_500))_COMPARISON_CELL_ENV_MARKED.txt',[1 1],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','JoinFiles/OtherStudies/Notepad Lists/TEMP.txt',[1 1],[1 0]);
FishLines('JoinFiles/Transcriptomics/Masuda2009_fin.txt','JoinFiles/Transcriptomics/TableS3.txt',[2 1],[1 1]);
FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/k12/(DB)TableS3_K12_ProteomeV13(for)TableS6_MatchingK12_BL21(DE3)_V2(1)NF.txt',[6 1],[1 0]);
FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2_k12once.txt','JoinFiles/OtherStudies/Notepad Lists/TEMP2.txt',[1 1],[1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%
WideStudies='JoinFiles/OtherStudies/Notepad Lists/InFiles_All studies.txt';
[CELLENV_wide_studies files_wide_studies CYTO_studies]=ExludeCyto(WideStudies,35,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/Protein List_Flotation_PRO008_1DG010_L1.txt',[],6,[2 2],[1 1]);

CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/CompareTrypticPeptides/Proteins_Large_Pho_Peptides_NonDet_CEP.txt',[],6,[2 1],[1 1]);

CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Matched.txt',[],5,[2 1],[1 1]);
CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Matched.txt',[],6,[2 4],[1 1]);
CellLoci('JoinFiles/k12/CatgID.TXT','JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Unique.txt',[],6,[2 2],[1 1]);


FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Common.txt',[6 2],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/mRNA_Taniguchi.txt',[2 1;3 1;4 2],[1 1]);
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/mRNA_present.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Matched.txt',[6 1],[1 1]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Matched.txt',[2 4],[1 1]);
FishLines('JoinFiles/BL2119/uniprot-organism_469008.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011CEPNonDet_BL21.txt',[6 4],[1 1]);

[list,Header,files,discr]=JoinMSMSResults('JoinFiles/k12/temp2.txt',1);
[list,Header,files,discr]=JoinMSMSResults('JoinFiles/k12/temp.txt',1);
FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2_k12once.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Common.txt',[1 2],[1 1]);

FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 2],[1 1]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Notfound.txt',[1 1],[1 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FishLines('JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 2],[1 1]);
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011.txt',[2 1],[1 1]);
FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011_Matched.txt',[6 1],[1 1]);
FishLines('JoinFiles/BL2119/Matching.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011.txt',[2 1],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FastaSubSelection('JoinFiles/k12/CellEnvk12.txt','Data/fASTA/ECOLI_K12_4407.fasta','acc',[1 0]);
FastaSubSelection('JoinFiles/BL2119/CellEnvBL21.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 0]);
%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/BL2119/Matching.txt','SignalP/BL21/NotFound.txt',[2 1],[1 1]);
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/uniprot-organism_83333.txt','SignalP/BL21/NotFound.txt',[6 1],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/k12/Ecoli_Sep2010_GO db.txt','SignalP/BL21/NotFound.txt',[2 1],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/k12/Ecoli_Sep2010_GO db.txt','JoinFiles/BL2119/Matching.txt',[3 1],[1 1]);
%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/BL2119/uniprot-organism_469008.txt','SignalP/BL21/BL21_lib.txt',[6 1],[1 1]);
%%%%%%%%%%%%%
PredictionToolsResults('SignalP/BL21/BL21_lib_acc.txt','SignalP/BL21/LipoShort_Unique.txt','SignalP/BL21/SignalBoth_Unique.txt','SignalP/BL21/TMHMM_Unique.txt','SignalP/BL21/Tat_Unique.txt','SignalP/BL21/Phobius_Unique.txt',1,1)
%%%%%%%%%%%%%%
FastaSubSelection('JoinFiles/BL2119/UniqueFV_acc.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 0]);
%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/BL2119/uniprot-organism_469008.txt','JoinFiles/BL2119/UniqueFV.txt',[6 1],[1 0]);
%%%%%
[list,Header,files,discr]=JoinMSMSResults('JoinFiles/Mpariam/Files.txt',1);
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/Mpariam/ColB.txt','JoinFiles/Mpariam/ColA.txt',[1 1],[1 1]);
%%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/DetectedTranscr.txt',[1 1],[1 1]);
CompareExpTheoPeptides('JoinFiles/TheoreticalPeptides/Theoretical.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/DetectedProtein.txt',[1 2 1 2 1 2],[1 1 1],1);
%%%%%%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/Oberto2009_.txt',[2 1;3 1;4 1],[1 1]);
%%%%%%%%%
WideStudies='JoinFiles/OtherStudies/Notepad Lists/Studies_wide_10.txt';
LB_All='JoinFiles/OtherStudies/Medium/LB_Merge.txt';
ExpressionAnalysis('JoinFiles/Transcriptomics/Merge_Protein_LB_temp.txt','JoinFiles/Transcriptomics/Merge_mRNA_LB_temp.txt',4);
ExpressionAnalysis('JoinFiles/Transcriptomics/Merge_Protein_LB.txt','JoinFiles/Transcriptomics/Merge_mRNA_LB.txt',3);
%%%%%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('Data/DataSetGO/DB_NoveV12.txt','JoinFiles/OtherStudies/Notepad Lists/OurStudyLoose.txt',[2 1],[1 0]);
%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('Data/DataSetGO/DB_NoveV12.txt','Data/TrypticPeptides/COUNTexp.txt',[2 1],[1 1]);
%%%%%%%%%%%
CompareExpTheoPeptides('JoinFiles/TheoreticalPeptides/Theoretical.txt','JoinFiles/TheoreticalPeptides/Experimental.txt','JoinFiles/Transcriptomics/NotDetByUs_mRNALevel(1)L.txt',[1 2 1 2 1 2],[1 1 1],1);
%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/mRNA_present.txt',[1 1],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/Wang2007.txt',[4 2],[1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%
Studies=[{'Molloy1999'}	{'Molloy2000'}	{'Gevaert2002'}	{'Yan2002'}	{'Fountoulakis2003'}	{'Corbin2003'}	{'Taoka2004'}	{'Butland2005'}	{'Lopez2005'}	{'Spelbrink2005'}	{'Stenberg2005'}	{'Wagner2005'}	{'Arifuzzaman2006'}	{'Baars2006'}	{'Huang2006'}	{'Ji2006'}	{'Lasserre2006'}	{'Marani2006'}	{'Cirulli2007'}	{'Wagner2007'}	{'Zhang2007'}	{'Jarchow2008'}	{'Wagner2008'}	{'Xia2008'}	{'Qian2008'}	{'Iwasaki2009'}	{'Masuda2009'}	{'Vertommen2009'}	{'Hemm2010'}	{'Iwasaki2010'}	{'Peng2010'}	{'Price2010'}	{'Wickstrom2010'}	{'OurStudy'}];
CEP=[5	36	311	64	146	403	480	491	232	46	53	159	733	158	59	394	107	6	63	103	353	34	146	64	21	840	800	64	20	840	42	71	52	778];
TotalCep=1935;
CEP_perc=(CEP./TotalCep)*100;
f=figure;[freq x]=hist(CEP_perc,0:5:100);bar(x,freq);
xlabel('Percent of CEP detected');
ylabel('# of Studies');
title('34 Studies including Us');

FileWriteTable('OtherStudies/Studies_wide_10.txt',Studies(CEP_perc>10)',[],'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%
Taniguchi=ReadTable('Transcriptomics/Taniguchi2010_fin.txt');
ppm=CellTable2Double(Taniguchi(2:end,end));
proteins=size(ppm,1);
MW=zeros(proteins,1);
for i=1:1:proteins
    [MW(i)]=MolMass(Taniguchi{i+1,2},'m');
end

molarity=( ppm )./MW;

FileWriteTable('Transcriptomics/Taniguchi2010_fin_mmol.txt',[[Taniguchi(1,:) {'mMol'}];[Taniguchi(2:end,:) Double2CellTable(molarity)]],[],'w');

%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','OtherStudies/Temp/ppm_Taniguchi2010.txt',[4 2],[1 1]);
%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','OtherStudies/reolderproteomicsstudy/Iwasaki2010.txt',[2 2;2 3;5 1],[1 1]);
%%% Match Peng table with transcritptomic and APEX data   %%%%
[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','Transcriptomics/PengTable.txt',[2 2;3 2;4 3],[1 0]);
%% THein Paper Merge
[list,Header]=JoinMSMSResults('OtherStudies/Temp/MergeMethods.txt',1);
%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','OtherStudies/Temp/Matched.txt',[6 1],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','OtherStudies/Temp/(DB)uniprot-organism_469008(for)TheinUniprot_NotFound(5)F.txt',[2 6],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('BL2119/uniprot-organism_469008.txt','OtherStudies/Temp/TheinUniprot_NotFound.txt',[5 1],[1 0]);
%% Fish Uniprot Acc Masuda 2009
% [FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','Transcriptomics/Masuda_et_al_Supple_rev2009.txt',[5 1;2 3;3 3],[1 1]);
% [FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','Transcriptomics/Ishihama2008_s2.txt',[1 1],[1 1]);
%% Fish GI for Thein Paper
TheoReadTable_=ReadTable('OtherStudies/Temp/ORList_MergeMethods.txt');
TheoReadTable=TheoReadTable_(2:end,:);
acc_expr='[A-Z][0-9A-Za-z]{5}(-[0-9])*';
locus_expr='[A-Z]+_[0-9]+';
locus_tag_expr=['locus_tag="',locus_expr];
exp=['Swiss-Prot:',acc_expr];
SWISSACC=cell(size(TheoReadTable,1),1);
FileWriteTable('OtherStudies/Temp/TheinUniprot.txt',TheoReadTable_(1,:),[],'w');
for i=607:1:size(TheoReadTable,1)
%     [GOT]=getgenpept(TheoReadTable{i,1});
    getgenpept(TheoReadTable{i,1}, 'ToFile','TempRetr.txt');
    GenPeptData =ReadTable('TempRetr.txt');
%     [CellTable]=CharTable2Cell(GOT.Features);
    [AccLine_st AccLine_en AccLine_ext AccLine_]=regexp(GenPeptData,locus_tag_expr); % Gene Accessions
    [AccLine]=CellTable2StrTable(AccLine_);
    [Acc_st Acc_en Acc_ext Acc]=regexp(AccLine,locus_expr); % Gene Accessions
    [LineFound]=find(CellTable2Double(Acc_st)>0);
    SWISSACC(i)=Acc{LineFound};
    FileWriteTable('TempRetr.txt',{'Hey'},[],'w');
    FileWriteTable('OtherStudies/Temp/TheinUniprot_NotFound.txt',[Acc{LineFound} TheoReadTable(i,:)],[],'a');
    pause(2);
end
%%%  Hyperlinks of Accession Numbers  %%%%
clear all;
TheoReadTable=ReadTable('k12/temp_All_K12.txt');
EmptyCell=cell(size(TheoReadTable,1),1);
% links{p}=['=HYPERLINK("http://www.uniprot.org/uniprot/?query=',IDsReadTable{p},'&sort=score")'];

Hyper_Links_nterm=MergeColumns([EmptyCell TheoReadTable],'=HYPERLINK("http://www.uniprot.org/uniprot/?query=');
Hyper_Links_mid=MergeColumns([Hyper_Links_nterm TheoReadTable],'&sort=score";"');
Hyper_Links_cterm=MergeColumns([Hyper_Links_mid EmptyCell],'")');
FileWriteTable('k12/Hyperlinks.txt',Hyper_Links_cterm,[],'w');
%%%%%%%%%%%%%%%%%%%%%%%%
[FOUNDTABLE count_found count_nfound]=FishLines('BL2119/uniprot-organism_469008.txt','BL2119/temp_allBL21.txt',[6 1],[1 0]);

[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','../Data/DataSetGO/K12_ProteomeV11.txt',[1 2],[1 1]);

[FOUNDTABLE count_found count_nfound]=FishLines('../Data/DataSetGO/K12_ProteomeV11.txt','k12/CoreProteome/2370protComm25str.txt',[2 1],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/24strCommon_withAcc_NotFound.txt',[2 1],[1 0]);
[FOUNDTABLE count_found count_nfound]=FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/24strCommon.txt',[31 1;32 1],[1 1]);
%%%
[FOUNDTABLE count_found count_nfound]=FishLines('../Data/DataSetGO/DB_NoveV10.txt','k12/CoreProteome/K12.txt',[2 1],[1 1]);
FileWriteTable('k12/CoreProteome/K12.txt',FOUNDTABLE,[],'w');
%%%%
FishLines('k12/CoreProteome/BL_Common.txt','k12/CoreProteome/CommonCore_MG1655_ref.txt',[18 1],[0 1]);
[FOUNDTABLE count_found count_nfound]=FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/CommonCore_MG1655_ref.txt',[1 1],[1 1]);
%%% Compare K12 remaining Common with BL %%%
[list,Header]=JoinMSMSResults('k12/CoreProteome/K12_VS_BL_Common.txt');
Table=CellTable2Double(list(:,2:end-1));
sum(Table(:,3)==1)
sum(Table(:,1)==1)
sum(Table(:,2)==1)

sum(Table(:,1)==1 & Table(:,2)==0)
sum(Table(:,2)==1  & Table(:,1)==0)
%% Prophage Common  %%%%%%
FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/BLCom_XOR_K12Com.txt',[1 1],[1 0]);
FishLines('k12/CoreProteome/K12_removed.txt','k12/CoreProteome/BLCom_XOR_K12Com.txt',[2 1],[0 0]);
[list,Header]=JoinMSMSResults('k12/CoreProteome/FishComonBLonly_from_K12removed.txt');
%%% Rename Insertion elements K12
FishLines('k12/K12Unique/K12Unique_NoveV10.txt','k12/InsertionElem/ORList_Files.txt',[2 1],[0 1]);
[list,Header,files]=JoinMSMSResults('k12/InsertionElem/Files.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% From b numbers to Uniprot Acc  %%%
% [FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','k12/CoreProteome/CommonCore_MG1655_ref.txt',[4 1],[1 1]);
% FileWriteTable('k12/CoreProteome/CommonCore_MG1655_ref.txt',FOUNDTABLE,[],'w');

%% Only in Eco Gene Core  %%%
FishLines('k12/K12Unique/K12Unique_NoveV9.txt','k12/CoreProteome/OnlyInEcoGeneCore.txt',[2 1],[1 0]);

%% Only in Eco Gene Core  %%%
[list,Header]=JoinMSMSResults('k12/InsertionElem/Files.txt');
%%%%
FishLines('k12/K12Unique/K12Unique_NoveV9.txt','k12/CoreProteome/CommonCore_MG1655_ref.txt',[4 1],[1 1]);
[FOUNDTABLE count_found count_nfound]=FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/CommonCore_MG1655_ref.txt',[1 1],[1 1]);
%%% Compare K12 remaining Common with BL %%%
[list,Header]=JoinMSMSResults('k12/CoreProteome/K12_VS_BL_Common.txt');
Table=CellTable2Double(list(:,2:end-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[list,Header]=JoinMSMSResults('k12/CoreProteome/CoreEcoGeneVSMultiGenome.txt');
Table=CellTable2Double(list(:,2:end-1));
sum(Table(:,3)==1)
sum(Table(:,1)==1)
sum(Table(:,2)==1)

sum(Table(:,1)==1 & Table(:,2)==0)
sum(Table(:,2)==1  & Table(:,1)==0)
list(Table(:,1)==1 & Table(:,2)==0,1)
%%  MultiGenomics K12 Unique 357  %%
FishLines('Data/DataSetGO/DB_SeptV5.txt','JoinFiles/K12Removed/Removed.txt',[2 1],[1 0]);

%% K12 Unique Fasta  %%
FastaSubSelection('JoinFiles/k12/K12Unique/K12Unique_NoveV9.txt','Data/fASTA/ECOLI_K12_4407.fasta','acc',[2 0]);

% Common Core %%%
[FOUNDTABLE count_found count_nfound]=FishLines('uniprot-organism_83333.txt','k12/CoreProteome/K12.txt',[1 1],[1 0]);

FileWriteTable('k12/CoreProteome/K12.txt',FOUNDTABLE(:,1:end-1),[],'w');

%% write file names
FILES=[{['k12/CoreProteome/CommonCore.txt',9,'4',9,'1']};{['k12/CoreProteome/K12.txt',9,'1',9,'1']}];
FileWriteTable('k12/CoreProteome/SelCommonCore.txt',FILES,[],'w');

% FILES=[{['k12/CoreProteome/CommonCore.txt',9,'5',9,'1']};{['k12/CoreProteome/CommonCore_MG1655_ref.txt',9,'1',9,'1']};{['k12/CoreProteome/K12.txt',9,'5',9,'1']}];
% FileWriteTable('k12/CoreProteome/CoreEcoGeneVSMultiGenome.txt',FILES,[],'w');

[list,Header]=JoinMSMSResults('k12/CoreProteome/SelCommonCore.txt');
Table=CellTable2Double(list(:,2:end-1));

[FOUNDTABLE count_found count_nfound]=FishLines('k12/CoreProteome/K12.txt','k12/CoreProteome/CommonCore.txt',[1 4],[1 1]);

sum(Table(:,1)==1 & Table(:,2)==1)
