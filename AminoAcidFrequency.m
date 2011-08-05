function [occur all_index]=AminoAcidFrequency(DataFile,len,step)

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(DataFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
AMINOS={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};

if(SPLen(1)==0)
    SPLen=SPLen+1;
end

SEQ=SEQ(TotLen>(len+SPLen));
Names=GnNames(TotLen>(len+SPLen));
SPLen=SPLen(TotLen>(len+SPLen));
Peptides=length(Names);
occur=zeros(Peptides,20);
% distr=zeros(20,len);

all_index=zeros(20,floor(len/step));

for j=1:1:20
    for i=1:1:Peptides
        [Pho start_idx end_idx matches]=CountPho(SEQ{i}(SPLen(i)+1:SPLen(i)+len),AMINOS{j});
%         cur_index=zeros(1,len);cur_index(start_idx)=1;
%         distr(j,:)=distr(j,:)+cur_index;
        occur(i,j)=Pho;
        [freq x]=hist(start_idx,step:step:len);
        all_index(j,:)=all_index(j,:)+freq;
    end
end

end