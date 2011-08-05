function [CellTable]=CharTable2Cell(CharTable)

[tokens columns]=size(CharTable);
CellTable=cell(tokens,1);

for i=1:1:tokens
    CellTable(i)={CharTable(i,:)};
end
end