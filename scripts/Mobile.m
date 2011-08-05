clear all;
close all;

[Table_Mobile]=ReadTable('JoinFiles/Mobile Elements/700Mobile.txt','\n');
[MobileElements]=CellTable2Double(Table_Mobile(2:end,2:end));
[Table_DB]=ReadTable('JoinFiles/Mobile Elements/EcoData022711-074325.txt','\n');
[Gene_DBS]=CellTable2Double(Table_DB(2:end,6:7));

Proteins=size(Gene_DBS,1);
index=zeros(Proteins,1);

FileWriteTable('JoinFiles/Mobile Elements/K12_Mobile.txt',[Table_DB(1,:) Table_Mobile(1,1)],[],'w');

for i=1:1:Proteins
    Mobile_L=MobileElements(:,1);
    Mobile_R=MobileElements(:,2);
    Gene_L=Gene_DBS(i,1);
    Gene_R=Gene_DBS(i,2);
    condition_within = ( Gene_L>=Mobile_L & Gene_L<Mobile_R ) & (Gene_R>Mobile_L & Gene_R<=Mobile_R);
    condition_Nterm = ( Gene_L<Mobile_L & Gene_L<Mobile_R ) & (Gene_R>Mobile_L & Gene_R<=Mobile_R);
    condition_Cterm = ( Gene_L>=Mobile_L & Gene_L<Mobile_R ) & (Gene_R>Mobile_L & Gene_R>Mobile_R);
    totally=condition_within | condition_Nterm | condition_Cterm;
    if(sum(totally)>0)
        FileWriteTable('JoinFiles/Mobile Elements/K12_Mobile.txt',[Table_DB(i+1,:) Table_Mobile(logical([0;totally]),:)],[],'a');
    end
end

FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/Mobile Elements/K12_Mobile.txt',[2 2;3 2;1 3;5 4],[1 1]);
FishLines('JoinFiles/Mobile Elements/(DB)uniprot-organism_83333_SubLocFixed(for)K12_Mobile(2)F.txt','JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 2],[1 1]);

JoinMSMSResults('JoinFiles/k12/temp3.txt',1);

FishLines('JoinFiles/k12/TableS6_MatchingK12_BL21(DE3)_V2.txt','JoinFiles/Mobile Elements/400Mobile_Elements.txt',[1 1],[1 1]);