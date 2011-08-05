function [CellTable]=StrMatrix2CellMatrix(strMatrix)

[row]=size(strMatrix,1);
CellTable=cell(row,1);

for i=1:1:row
    CellTable(i)={strMatrix(i,:)};
end
end