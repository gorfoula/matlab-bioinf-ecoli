function [SEQ]=RandomSeq(peptides,len,outSubFolder)

FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};

index=round(rand(peptides,len)*19);  % random numbers from 0-19
index=index+1; % shift to make it index
index_=reshape(index,1,peptides*len);
SEQ_=FEATURES(index_);

SEQ=reshape(SEQ_,peptides,len);
[Comb_SEQ]=MergeColumns(SEQ);

FileWriteTable([outSubFolder,'RandomSeq.txt'],Comb_SEQ,'Sequence','w');
end