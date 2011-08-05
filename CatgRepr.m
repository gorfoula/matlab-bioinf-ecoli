function [Names,Values]=CatgRepr(Sequence,type,fileCat)

[Names,CatgCors,CATG]=LoadAAFeaturesCatg(fileCat);

[Row,Col,Zax]=size(Sequence);

switch (type)
     case 'int'
         Values=zeros(Row,Col,Zax);
         seqLen=Col;
     case 'bin'
         Values=zeros(Row,CATG,Zax);
         seqLen=Zax;
    otherwise
         display('Not rigth input for TYPE (int/char/bin)');
         return;
end

for peptide=1:1:Row
    for aa=1:1:seqLen
        switch (type)
                case 'int'
                    Values(peptide,aa)=CatgCors(Sequence(peptide,aa));
                case 'bin'
                    index=1:1:Col;
                    amino_type=index*Sequence(peptide,:,aa)';
                    if(amino_type>0)
                        Values(peptide,CatgCors(amino_type),aa)=1;
                    end
            otherwise
                    display('Not rigth input for TYPE (int/char/bin)');
                    return;
        end
    end
end

end