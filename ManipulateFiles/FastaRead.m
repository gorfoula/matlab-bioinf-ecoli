function [HEAD SEQ]=FastaRead(IFasta)

[lines] = textread(IFasta,'%s',-1,'delimiter','\n','bufsize',21000);
headers=CountHeaders(lines);
HEAD=cell(1,headers);
SEQ=cell(1,headers);
count=0;
for i=1:1:length(lines)
    [header_flag]=IsItHeader(lines{i});
    if(header_flag)
        HEAD{count+1}=lines{i};count=count+1;
    else
        SEQ{count}=[SEQ{count},lines{i}];
    end
end

end

function [count]=CountHeaders(lines)
count=0;
    for i=1:1:length(lines)
        found=regexp(lines{i},'>');
        count=count+IsItHeader(lines{i});
    end
end

function [header_flag]=IsItHeader(line)
    found=regexp(line,'>');
    header_flag=~isempty(found);
end