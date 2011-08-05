function []=SavefileSplitSPMature(IFile)

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

samples=length(GnNames);
match=regexp(IFile,'.txt|/');

s=fopen(['Data/',IFile(match(1)+1:match(2)-1),'split_SP_MAT.txt'],'w');
fprintf(s,'Names \t Description \t Total Len \t MW \t SP Len \t Sequence \t SP Sequence \t MATURE Sequence\n');
for ig=1:1:samples
        fprintf(s,'%s \t %s \t %d \t %d \t %d \t %s \t %s \t %s \n',GnNames{ig},Descr{ig},TotLen(ig),MW(ig),SPLen(ig),SEQ{ig},SEQ{ig}(1:SPLen(ig)),SEQ{ig}(SPLen(ig)+1:end));
end

fclose all;

end