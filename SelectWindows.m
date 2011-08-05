function []=SelectWindows(IFs,len)

[FileNames] = textread(IFs,'%s',-1,'delimiter','\t');

[GnNamesA, DescrA, TotLenA, MWA, SPLenA, SEQA] = textread(FileNames{1},'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNamesB, DescrB, TotLenB, MWB, SPLenB, SEQB] = textread(FileNames{2},'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

% if(SPLenA(1)>0)  % if the protein is a Secreted one we exclude Meth therefore the SP length is minus one
%     SPLen=SPLen-1;
% end

SELECTEDA=SEQA(TotLenA>=len+1);    %Select Petides with Seq long enough
samples=length(SELECTEDA);
%%%%% INIT %%%%%%
AAs_intA=zeros(samples,len);
SelSeqA=cell(samples,1);
PhoMatrixA=zeros(samples,len);
PhoSumA=zeros(samples,1);
BulkMatrixA=zeros(samples,len);
BulkSumA=zeros(samples,1);
PolMatrixA=zeros(samples,len);
PolSumA=zeros(samples,1);
PolSumMatrix_negA=zeros(samples,len);
PolSum_negA=zeros(samples,1);
for i=1:1:samples
        [AAs_intA(i,1:len) SelSeqA(i)]=AARepresentation('int',2,len+1,{SELECTEDA{i}});
        [PhoMatrixA(i,1:len) PhoSumA(i)]=Hydrophobicity(AAs_intA(i,1:len),'En');
        [BulkMatrixA(i,1:len) BulkSumA(i)]=Bulckiness(AAs_intA(i,1:len));
        [PolMatrixA(i,1:len) PolSumA(i)]=Polarity(AAs_intA(i,1:len),'Pos');
        [PolSumMatrix_negA(i,1:len) PolSum_negA(i)]=Polarity(AAs_intA(i,1:len),'Neg');
end
SELECTEDB=SEQB(TotLenB>=len+1);    %Select Petides with Seq long enough
samples=length(SELECTEDB);
%%%%% INIT %%%%%%
AAs_intB=zeros(samples,len);
SelSeqB=cell(samples,1);
PhoMatrixB=zeros(samples,len);
PhoSumB=zeros(samples,1);
BulkMatrixB=zeros(samples,len);
BulkSumB=zeros(samples,1);
PolMatrixB=zeros(samples,len);
PolSumB=zeros(samples,1);
PolSumMatrix_negB=zeros(samples,len);
PolSum_negB=zeros(samples,1);
for i=1:1:samples
        [AAs_intB(i,1:len) SelSeqB(i)]=AARepresentation('int',2,len+1,{SELECTEDB{i}});
        [PhoMatrixB(i,1:len) PhoSumB(i)]=Hydrophobicity(AAs_intB(i,1:len),'En');
        [BulkMatrixB(i,1:len) BulkSumB(i)]=Bulckiness(AAs_intB(i,1:len));
        [PolMatrixB(i,1:len) PolSumB(i)]=Polarity(AAs_intB(i,1:len),'Pos');
        [PolSumMatrix_negB(i,1:len) PolSum_negB(i)]=Polarity(AAs_intB(i,1:len),'Neg');
end


[MaxWindow_Pho MaxDiff_Pho]=AllWindowPropertySums(PhoMatrixA,PhoMatrixB,[2:1:70])
[MaxWindow_Bulk MaxDiff_Bulk]=AllWindowPropertySums(BulkMatrixA,BulkMatrixB,[2:1:70])
[MaxWindow_Pos MaxDiff_Pos]=AllWindowPropertySums(PolMatrixA,PolMatrixB,[2:1:70])
[MaxWindow_Pol MaxDiff_Pol]=AllWindowPropertySums(PolSumMatrix_negA,PolSumMatrix_negB,[2:1:70])

end
