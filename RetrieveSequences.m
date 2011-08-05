function []=RetrieveSequences(IFile,tableF,splits)

[EG GN SP GI JW ASAP EB Len Description MW]=textread(tableF,'%s %s %s %s %s %s %s %s %s %s',-1,'delimiter','\t');
[SwissProt] = textread(IFile,'%s',-1,'delimiter','\t');

Peptides=length(SwissProt);
index=zeros(Peptides,1);
s1 = regexp(IFile,'[.]');
for j=1:1:Peptides
    idx=find(strcmp(upper(GN),upper(SwissProt{j}))==1);
    if(length(idx)>0)
        index(j)=idx;
    else
        s=fopen([IFile(1:s1(1)-1),'Notfound.txt'],'a');
        fprintf(s,'%s\n',SwissProt{j});
        fclose(s);
    end
end


s=fopen([IFile(1:s1(1)-1),'.fasta'],'a');
for i=1:Peptides
    i
    if(index(i)>0)
        if(strcmp(SP{index(i)},'Null')==0)
            entry_name=upper(GN{index(i)});
            STRUCT = getgenpept(SP{index(i)});
            SEQ=upper(STRUCT.Sequence);
            [SLSQ]=SequenceSplit(SEQ,splits);
            fprintf(s,'%s\n',['>',entry_name,'_ECOLI|',SwissProt{i},'|',Description{i}]);
            for i=1:length(SLSQ)
                fprintf(s,'%s\n',SLSQ{i});
            end
            fprintf(s,'\n');
        else
            Sn=fopen('Dataset/Null.txt','a');
            fprintf(Sn,'%s\n',GN{index(i)});
            fclose(Sn);            
        end
        
    end
    
end

fclose(s);