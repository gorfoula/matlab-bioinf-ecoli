
ProteomeDB_=ReadTable('JoinFiles/BL2119/CommonCoreBL21_K12_V2.txt');
SubLoci=ProteomeDB_(2:end,20);
index_cyto=(strcmp(SubLoci,'A')) | (strcmp(SubLoci,'A_trl'));
index_celenv=not(index_cyto);

FileWriteTable('JoinFiles/BL2119/CellEnvelopeBL21_K12_V2.txt',ProteomeDB_(index_celenv,:),[],'w');

FastaSubSelection('JoinFiles/BL2119/CellEnvelopeBL21_K12.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 0])

FastaSubSelection('JoinFiles/BL2119/BL2119NovV2_NoMobileElement.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 1])
FastaSubSelection('JoinFiles/BL2119/CommonCoreBL21_K12_V2.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 1])
FastaSubSelection('JoinFiles/BL2119/UniqueBL21_V2.txt','JoinFiles/BL2119/uniprot-organism_469008.fasta','acc',[1 0])