function [Table index]=CellTable2Double(CellTable)
[tokens cols]=size(CellTable);
Table=zeros(tokens,cols);
index=ones(tokens,cols);
for i=1:1:tokens
    for j=1:1:cols
        if(isempty(CellTable{i,j})==0)
            if(isnumeric(CellTable{i,j}))
                Table(i,j)=CellTable{i,j};
            else
                Table(i,j)=str2double(CellTable{i,j});
                index(i,j)=not(isnan(Table(i,j)));
            end
        else
            index(i,j)=0;
        end
    end
end
index=logical(index);
end