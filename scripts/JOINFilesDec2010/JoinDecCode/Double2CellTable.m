function [CellTable]=Double2CellTable(DoubleTable)

[tokens columns]=size(DoubleTable);
if(columns==0)
    CellTable=cell(tokens,1);columns=1;
else
    CellTable=cell(tokens,columns);
end


for i=1:1:tokens
    for n=1:columns
        CellTable(i,n)={DoubleTable(i,n)};
    end
%     CellTable{i}=num2str(DoubleTable(i));
end
end