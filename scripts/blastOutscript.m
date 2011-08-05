clear all;clc;
[BlastOut] = ReadTable('JoinFiles/k12/K12Unique/BlastUnique_K12NoveV9.txt');

query_hit='[0-9]+[.]*[0-9]+';
NoHits_expr='#\s+[0-9]+\s+';
mnem_expr='[A-Za-z0-9]{3,6}_ECOLI';
embl_expr='[A-Z][A-Za-z0-9]+';
lines=size(BlastOut(:,1),1);

%%%%% ROWS WITH HITS %%%%%
[start_idx]=regexp(BlastOut(:,end),query_hit);
rows_hit=CellTable2Double(start_idx)>0;
%%%%% MATCH _ECOLI HITS %%%%%
[EcoHits_st EcoHits_en EcoHits_ext EcoHits_]=regexp(BlastOut(rows_hit,2),mnem_expr);
EcoHits=CellTable2StrTable(EcoHits_);
%%%%% MATCH EMBL IDS FROM HITS %%%%%
[EmblIDs_st EmblIDs_en EmblIDs_ext EmblIDs_]=regexp(BlastOut(rows_hit,2),embl_expr);
EmblIDs=CellTable2StrTable(EmblIDs_);
%%%% # of HITS OF EACH QUERY %%%%
[start_idx]=regexp(BlastOut(:,1),NoHits_expr);
rows=CellTable2Double(start_idx)>0;
[start_idx, end_idx, extents, matches]=regexp(BlastOut(rows,1),'[0-9]+');
No_Hits=CellTable2Double(matches);

BlastOut=BlastOut(rows_hit,:);
%%%%% MATCH _ECOLI HITS %%%%%
[Acc_st Acc_en Acc_ext Acc_]=regexp(BlastOut(:,1),'[0-9A-Z]{6}');
Acc=CellTable2StrTable(Acc_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
queries=size(No_Hits,1);
Identidy=CellTable2Double(BlastOut(:,3));
evalue=CellTable2Double(BlastOut(:,end-1));
product=Identidy.* ( 1./evalue);
count=1;i=1;
FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_Unique.txt',{'Hey'},[],'w');
FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_Matched.txt',{'Hey'},[],'w');
FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_AllWithHits.txt',{'Hey'},[],'w');
FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut.txt',[Acc EcoHits BlastOut],[],'w');
for i=1:1:queries
    count_end=count+No_Hits(i)-1;
    [count count_end]
    cur_products=product(count:count_end,:);
    [m pos]=max(cur_products);
    if(evalue(count+pos-1)>10^(-3))
        FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_Unique.txt',[BlastOut(count+pos-1,1:2) Acc(count+pos-1,1) EcoHits(count+pos-1,:) Identidy(count+pos-1) evalue(count+pos-1)],[],'a');
    else
        FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_Matched.txt',[BlastOut(count+pos-1,1:2) Acc(count+pos-1,1) EcoHits(count+pos-1,:) Identidy(count+pos-1) evalue(count+pos-1)],[],'a');
    end
    FileWriteTable('JoinFiles/BL2119/Unique blast/BlastOut_AllWithHits.txt',[BlastOut(count+pos-1,1:2) Acc(count+pos-1,1) EcoHits(count+pos-1,:) Identidy(count+pos-1) evalue(count+pos-1)],[],'a');
    count=count_end+1;    
end

FishLines('JoinFiles/BL2119/Unique blast/BlastOut_AllWithHits.txt','JoinFiles/BL2119/UniqueBL21_V2.txt',[1 1],[1 0]);
FishLines('JoinFiles/BL2119/Unique blast/BlastOut_Matched.txt','JoinFiles/BL2119/UniqueBL21_V2.txt',[1 1],[1 0]);
FishLines('JoinFiles/BL2119/Unique blast/BlastOut_Unique.txt','JoinFiles/BL2119/UniqueBL21_V2.txt',[1 1],[1 0]);
FishLines('JoinFiles/uniprot-organism_83333.txt','JoinFiles/BL2119/Unique blast/ComnK12.txt',[1 1],[1 0]);

FishLines('Data/DataSetGO/DB_NoveV8(NoInsElemPhageEcoRhs).txt','JoinFiles/BL2119/Unique blast/temp.txt',[1 1],[1 0]);
FishLines('Data/DataSetGO/DB_SeptV5.txt','JoinFiles/BL2119/Unique blast/(DB)DB_NoveV8(NoInsElemPhageEcoRhs)(for)temp(1)NF.txt',[1 1],[1 0]);
FishLines('Data/DataSetGO/DB_NoveV8(NoInsElemPhageEcoRhs).txt','JoinFiles/BL2119/Unique blast/temp.txt',[1 1],[1 0]);
% CellLoci('Data/DataSetGO/CatgID.TXT','Data/DataSetGO/DB_NoveV8(NoInsElemPhageEcoRhs).txt','JoinFiles/BL2119/Unique blast/temp.txt',[],6,[1 1],[1 1]);
