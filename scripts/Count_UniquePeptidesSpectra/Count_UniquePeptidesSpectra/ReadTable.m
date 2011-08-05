function [Table]=ReadTable(DatabaseIF,del)

if(exist('del','var')==0)      %if expr is not defined
    del='\n';
end

display(['Trying to read: ',DatabaseIF]);
[text] = textread(DatabaseIF,'%s',-1,'delimiter',del,'bufsize',10000);
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text(1:end),'[\t]');
Table=CellTable2StrTable(Table_);

end