function [Comb_SEQ]=MergeColumns(SEQ,separator,tailORhead)

if(exist('separator','var')==0)
   separator=9 ;
end
if(isempty(SEQ))
    Comb_SEQ=SEQ;
    return;
end
[peptides cols]=size(SEQ);
Comb_SEQ=SEQ(:,1);
if(cols==1)
    
    if (exist('tailORhead','var')==0)
%         display('Your input is a unicoloun Array. PLease specify were to attache the separator head|tail => h|t');
        tailORhead='t';
    end
    switch (tailORhead)
        case 't'
            Comb_SEQ=cellfun(@(a) [a,separator],SEQ,'uni',false);
        case 'h'
            Comb_SEQ=cellfun(@(a) [separator,a],SEQ,'uni',false);
        otherwise
            display('MergeColoumns line 19: tailORhead parameter should be t|h');
            return;
    end
else
    for i=2:cols
        Comb_SEQ=cellfun(@(a,b) [a,separator,b],Comb_SEQ,SEQ(:,i),'uni',false);
    end
end

end