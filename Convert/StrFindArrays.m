function [COUNT EXACT ALL]=StrFindArrays(ARRAY,str)

rows=size(ARRAY,1);
EXACT=zeros(rows,1);
ALL=zeros(rows,1);
COUNT=0;
for i=1:rows
    [StPos EnPos]=regexp(ARRAY{i},str);
    if(isempty(StPos)==0)
        ALL(i)=1;
        COUNT=COUNT+1;
        sizeDiff=size(ARRAY{i},2)-size(str,2);
        if(sizeDiff==0)
            EXACT(i)=1;
        end
    end
end

end