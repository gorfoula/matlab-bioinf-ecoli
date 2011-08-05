function []=FastaWrite(IFasta,HEAD,SEQ,split,open)
display(['fasta dir: ',IFasta]);
if(exist('open','var')==0)
    open='w';
end
s=fopen(IFasta,open);  %% file save group of proteins with specific best SP
headers=length(HEAD);
count=0;
for i=1:1:headers
    if(isempty(HEAD{i})==0 && isempty(HEAD{i})==0)
        if(strcmp(HEAD{i}(1),'>'))
            fprintf(s,[HEAD{i},'\n']);
        else
            fprintf(s,['>',HEAD{i},'\n']);
        end
        PrintSequence(s,SEQ{i},split);
    end
end
fclose(s);
end

function []=PrintSequence(s,SEQ,split)
len=length(SEQ);
lines=floor(len/split);
for i=0:1:lines-1
    fprintf(s,[SEQ(i*split+1:(i+1)*split),'\n']);
end
fprintf(s,[SEQ(lines*split+1:end),'\n\n']);
end