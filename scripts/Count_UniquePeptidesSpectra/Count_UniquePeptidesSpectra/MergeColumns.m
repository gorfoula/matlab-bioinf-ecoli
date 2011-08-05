function [Comb_SEQ]=MergeColumns(SEQ,separator)

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
    Comb_SEQ=cellfun(@(a) [a,separator],SEQ,'uni',false);
else
    for i=2:cols
        Comb_SEQ=cellfun(@(a,b) [a,separator,b],Comb_SEQ,SEQ(:,i),'uni',false);
    end
end

end