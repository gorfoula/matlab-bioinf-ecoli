% close all;clc;clear all;

VenDiagram(RandLevels(10,2));

[list,Header]=JoinMSMSResults('K12Unique/Files.txt');

mnem_expr='[A-Za-z0-9]{3,6}_ECOLI';
TheoReadTable=ReadTable('JoinFiles/BL2119/Unique blast/Mapped_EMBL.txt');
[start_idx , end_idx, extents, matches_]=regexp(TheoReadTable(:,3),mnem_expr);
matches=CellTable2StrTable(matches_);
rows=CellTable2Double(start_idx)>0;
FileWriteTable('JoinFiles/BL2119/Unique blast/Mapped_EMBL_ECOLI.txt',[matches(rows,:) TheoReadTable(rows,:)],[],'w');


FishLines('Data/DataSetGO/DB_SeptV5.txt','JoinFiles/K12Removed/Removed.txt',[2 1],[1 0]);
FishLines('JoinFiles/K12Removed/Removed.txt','Data/DataSetGO/DB_NoveV8(NoInsElemPhageEcoRhs).txt',[1 2],[0 1]);


FishLines('JoinFiles/K12Unique/Common.txt','JoinFiles/K12Unique/K12Total.txt',[1 2],[0 0]);

FishLines('JoinFiles/K12Unique/BL21Remv.txt','JoinFiles/K12Unique/Douplicates.txt',[17 1],[0 0]);

FastaSubSelection('JoinFiles/K12Unique/DouplicatesBL.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 0]);

FishLines('JoinFiles/BL2119/Unique blast/BlastOut.txt','JoinFiles/BL2119/Unique blast/Mapped_EMBL_ECOLI.txt',[1 1],[0 1]);

FishLines('JoinFiles/BL2119/Unique blast/BlastOut.txt','JoinFiles/BL2119/Unique blast/DouplicatesBL.txt',[1 1],[0 0]);

FishLines('JoinFiles/K12Unique/BL21Remv.txt','JoinFiles/K12Unique/Unique blast/temp.txt',[1 1],[0 0]);
