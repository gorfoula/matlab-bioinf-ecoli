function [Table]=CellTable2Double(CellTable)
[tokens cols]=size(CellTable);
Table=zeros(tokens,cols);

for i=1:1:tokens
    for j=1:1:cols
        if(isempty(CellTable{i,j})==0)
            if(isnumeric(CellTable{i,j}))
                Table(i,j)=CellTable{i,j};
            else
                Table(i,j)=str2double(CellTable{i,j});
            end
        end
    end
end