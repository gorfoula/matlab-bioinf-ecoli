function [count]=PrintCell(Ofile,CellTable,write_mode)
file_indx=fopen(Ofile,write_mode);
[proteins coloumns]=size(CellTable)
count=0;
for p=1:1:proteins
    for c=1:1:coloumns
        fprintf(file_indx,'%s\t',CellTable{p,c});  %% file save best SP in an index
    end
    fprintf(file_indx,'\n');
    count=count+1;
end

end