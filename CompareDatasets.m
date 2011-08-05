function []=CompareDatasets(File1,File2)

found=regexp(File1,'[/]');name=regexp(File1,'[.]');
dirdataset=File1(1:found(end));
fileName1=File1(found(end)+1:name(end)-1);
found=regexp(File2,'[/]');name=regexp(File2,'[.]');
fileName2=File2(found(end)+1:name(end)-1);

[GnNames1, Descr1, TotLen1, MW1, SPLen1, SEQ1] = textread(File1,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNames2, Descr2, TotLen2, MW2, SPLen2, SEQ2] = textread(File2,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
len=40;
sub_SEQ1=SubStrings(SEQ1,len);
sub_SEQ2=SubStrings(SEQ2,len);
for i=1:1:length(GnNames2)
    indx_seq_sub=find(strcmpi(sub_SEQ1,sub_SEQ2{i}));
    indx_total=find(strcmpi(SEQ1,SEQ2{i}));
    indx_name=find(strcmpi(GnNames1,GnNames2{i}));
    if(isempty(indx_name))
        AppendDataset([dirdataset,'NameDiff(',fileName1,'-',fileName2,').txt'],GnNames2(i), Descr2(i), TotLen2(i), MW2(i), SPLen2(i), SEQ2(i),'a');
    end
    if(length(indx_total)>1)
        AppendDataset([dirdataset,'DoublesTotalSeq(',fileName1,'-',fileName2,').txt'],GnNames2(i), Descr2(i), TotLen2(i), MW2(i), SPLen2(i), SEQ2(i),'a');
    elseif(length(indx_name)>1)
        AppendDataset([dirdataset,'DoublesNames(',fileName1,'-',fileName2,').txt'],GnNames2(i), Descr2(i), TotLen2(i), MW2(i), SPLen2(i), SEQ2(i),'a');
    elseif(length(indx_seq_sub)>1)
        AppendDataset([dirdataset,'DoublesSubSeq(',fileName1,'-',fileName2,').txt'],GnNames2(i), Descr2(i), TotLen2(i), MW2(i), SPLen2(i), SEQ2(i),'a');
    end
    if(isempty(indx_name))
        AppendDataset([dirdataset,'NotFound(',fileName1,'-',fileName2,').txt'],GnNames2(i), Descr2(i), TotLen2(i), MW2(i), SPLen2(i), SEQ2(i),'a');
    end
end
end

function [SubStrings]=SubStrings(StringList,len)
STRINGS=length(StringList);
SubStrings=cell(STRINGS,1);
for i=1:1:STRINGS
    if(length(StringList{i})<len)
        SubStrings{i}=StringList{i}(1:end);
    else
        SubStrings{i}=StringList{i}(1:len);
    end
end

end

function []=AppendDataset(oFile,GnNames, Descr, TotLen, MW, SPLen, SEQ,mode)
s=fopen(oFile,mode);  %% file save group of proteins with specific best SP
for i=1:1:length(GnNames)
    if(isempty(SEQ{i})==0)
        fprintf(s,[GnNames{i},'\t',Descr{i},'\t',num2str(TotLen(i)),'\t',num2str(MW(i)),'\t',num2str(SPLen(i)),'\t',SEQ{i},'\n']);
    end
end
fclose all;
end

