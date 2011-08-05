function [CellTable]=Double2CellTable(DoubleTable,sel)

if(exist('sel','var')==0)
    sel=0;
end

[tokens columns]=size(DoubleTable);
if(columns==0)
    CellTable=cell(tokens,1);columns=1;
else
    CellTable=cell(tokens,columns);
end


for i=1:1:tokens
    for n=1:columns
        if(sel==1)
            CellTable(i,n)={num2str(DoubleTable(i,n))};
        else
            CellTable(i,n)={DoubleTable(i,n)};
        end
        
    end
%     CellTable{i}=num2str(DoubleTable(i));
end
end