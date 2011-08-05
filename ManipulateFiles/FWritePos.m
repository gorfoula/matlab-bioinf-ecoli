function []=FWritePos(fid,Table)
if(isempty(Table)==1)
   fclose all;
   return;
end
peptides=length(Table(:,1));
column=length(Table(1,:));
for j=1:1:peptides   %for all best signal peptides
    for c=1:1:column-1
        fprintf(fid,'%s\t',num2str(Table{j,c}));  %% file save best SP in an index
    end
    fprintf(fid,'%s\n',num2str(Table{j,column}));  %% file save best SP in an index
end
end