function [Table index]=CellTable2StrTable(CellTable)
tokens=length(CellTable);
[rows columns]=size(CellTable{1,1});
if(columns==0)
    Table=cell(tokens,1);columns=1;
else
    Table=cell(tokens,columns);
end
index=zeros(tokens,1);
for i=1:1:tokens
   if(isempty(CellTable{i})==0)
       c=size(CellTable{i},2);
       if(c<columns)
           Table(i,1:c)=CellTable{i}(1:c);
       elseif(c>columns)
           Table=[Table cell(tokens,c-columns)];
           columns=c;
           Table(i,1:columns)=CellTable{i}(1:columns);
       else
           Table(i,1:columns)=CellTable{i}(1:columns);
       end
       index(i)=1;
   end
end

index=find(logical(index));
end