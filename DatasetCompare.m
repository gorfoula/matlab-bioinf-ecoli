function []=DatasetCompare(file1,file2)

found=regexp(file1,'[/]');name=regexp(file1,'[.]');
filename1=file1(found(end)+1:name(end)-1);
found2=regexp(file2,'[/]');name2=regexp(file2,'[.]');
dirdataset=file2(1:found2(end));
filename2=file2(found2(end)+1:name2(end)-1);

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(file1,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNames2, Descr2, TotLen2, MW2, SPLen2, SEQ2] = textread(file2,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

peptides=length(GnNames2);

for i=1:1:peptides
    indx_seq=find(strcmpi(SEQ,SEQ2{i}));
    if(isempty(indx_seq)==0)
        AppendDataset([dirdataset,'CommonProteins(',filename1,'VS',filename2,').txt'],{GnNames{indx_seq}}, {Descr{indx_seq}}, TotLen(indx_seq), MW(indx_seq), SPLen(indx_seq), {SEQ{indx_seq}},'a');
    end
end

end

function []=AppendDataset(oFile,GnNames, Descr, TotLen, MW, SPLen, SEQ,mode)
s=fopen(oFile,mode);  %% file save group of proteins with specific best SP
for i=1:1:length(GnNames)
    fprintf(s,[GnNames{i},'\t',Descr{i},'\t',num2str(TotLen(i)),'\t',num2str(MW(i)),'\t',num2str(SPLen(i)),'\t',SEQ{i},'\n']);
end
fclose(s);
end