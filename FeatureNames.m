function [FEATURES]=FeatureNames(ftsel,DomDimen)

switch (ftsel)
    case 20
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};
    case 103
        FEATURES=cell(1,sum(DomDimen));
        index=conv(DomDimen,ones(1,length(DomDimen)));
        FEATURES(1:index(1))=Concatenate('P',DomDimen(1));
        FEATURES(index(1)+1:index(2))=Concatenate('H',DomDimen(2));
        FEATURES(index(2)+1:index(3))=Concatenate('Ppos',DomDimen(3));
        FEATURES(index(3)+1:index(4))=Concatenate('Hpos',DomDimen(4));
        
    otherwise
        display('=>Wrong feature selection type');
        return;
end

end

function [fnames]=Concatenate(letter,fnum)

fnames=cell(1,fnum);

for i=1:1:fnum
    fnames{i}=[letter,num2str(i)];
end

end