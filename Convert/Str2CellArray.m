function [CellTable]=Str2CellArray(str)

[tokens]=length(str);
CellTable=cell(tokens,1);

for i=1:1:tokens
    CellTable(i)={str(i)};
end
end