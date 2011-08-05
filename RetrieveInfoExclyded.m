function []=RetrieveInfoExclyded(IFile)

[GeneNames] = textread(['Gems Results/Excluded/',IFile],'%s',-1,'delimiter','\t');

[eck bnum genbank gene	type product] = textread('geneproductfunctions.txt','%s %s %s %s %s %s',-1,'delimiter','\t');

Peptides=length(GeneNames)
Ecoli_Peptides=length(gene);

index=ones(1,Peptides);
not_index=zeros(1,Peptides);

for i=1:1:Peptides
    for j=1:1:Ecoli_Peptides
        if(strcmp(gene{j},GeneNames{i})==1 && strcmp(genbank{j},'')==0)
            index(i)=j;
            flag=1;
        end
    end
    if(flag==0)
        not_index(i)=i;
    end
    flag=0;
end

index=index(index>1);
not_index=not_index(not_index>0);
not_GeneNames=GeneNames(not_index);

GN=gene(index);
GeneIDs=genbank(index);
found=length(index)
Sequence=cell(found,1);
Definition=cell(found,1);
Length=zeros(found,1);

for i=1:1:length(index)
    s=fopen(['Gems Results/Excluded/',IFile,'.test.txt'],'a');
    GN{i}
    Seq = getgenpept(GeneIDs{i});
    getgenpept(GeneIDs{i},'ToFile',['Gems Results/Excluded/',IFile,'.fasta'],'FileFormat','FASTA');
    Sequence{i}=Seq.Sequence;
    Definition{i}=Seq.Definition;
    Length(i)=str2double(Seq.LocusSequenceLength);
    GeneMW=molweight(Sequence{i});
    
    fprintf(s,'%s\t%s\t%s\t%d\t%d\t%d\t%s\n',GeneIDs{i},GN{i},Definition{i},Length(i),GeneMW,1,upper(Sequence{i}));
    fclose(s);
    pause(5);
end
s=fopen(['Gems Results/Excluded/',IFile,'_xxxxx.txt'],'a');
for i=1:1:length(not_index)
    fprintf(s,'%s\n',not_GeneNames{i});
end
fclose all;
end


% IMP cyclohydrolase
% IMP synthetase
% IMP--aspartate ligase
% IMPase