function [Domains fnames max_Nend max_Hend]=DomainsLoad(Names,IFile)

[GeneNames Nend Hend Cend Nlen Hlen Clen] = textread(IFile,'%s %d %d %d %d %d %d',-1,'delimiter','\t');

max_Nend=max(Nend);
max_Hend=max(Hend);

Peptides=length(Names);
Domains=zeros(Peptides,6);
fnames=cell(Peptides,1);

for i=1:1:Peptides
%     if(i==406)
%         display('Hey baby');
%     end
    indx=find(strcmpi(GeneNames,Names{i}));
    if(isempty(indx)==0)
        Domains(i,:)=[Nend(indx) Hend(indx) Cend(indx) Nlen(indx) Hlen(indx) Clen(indx)];
        fnames{i}=Names{i};
    else
        display([Names{i},' Not found']);
    end
end

end