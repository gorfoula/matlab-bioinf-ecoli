function [positives Pho PhoSum pos_indx pho_indx DomSeq]=CountPropertiesDomain(SQ,SQNames,Names,Domains,max_Nend,max_Hend)

Peptides=length(Names);

positives=ones(Peptides,1)*(-100);pos_indx=zeros(Peptides,max_Nend);
Pho=ones(Peptides,1)*(-100);pho_indx=zeros(Peptides,max_Hend);

PhoSum=ones(Peptides,1)*(-100);
PhoMatrix=zeros(Peptides,12);
DomSeq=cell(Peptides,3);
for i=1:1:Peptides
    ind=strcmp(SQNames,Names{i});
    if(sum(ind)>0)
        if(length(SQ{ind})<Domains(i,3))
            display([Names{i},' Larger than SP']);
        else
            NDomainSeq=SQ{ind}(1:Domains(i,1));
            HDomainSeq=SQ{ind}(Domains(i,1)+1:Domains(i,2));
            CDomainSeq=SQ{ind}(Domains(i,2)+1:end);
            DomSeq(i,:)={NDomainSeq HDomainSeq CDomainSeq};
            Hlen=length(HDomainSeq);
            [positives(i) indx]=CountPho(NDomainSeq,'[KR]');pos_indx(i,indx)=1;
            [Pho(i) indx]=CountPho(HDomainSeq);pho_indx(i,indx+Domains(i,4))=1;

            [AAs_int temp]=AARepresentation('int',1,Hlen,{HDomainSeq});
            [PhoMatrix(i,1:Hlen) PhoSum(i)]=Hydrophobicity(AAs_int,'En');
        end
    else
        display([Names{i},' Not found']);
    end
end

end