%% Write file with proteins assigned to categories
function []=FileWriteTable(OFile,Table,header,write_mode)
[file_indx]=fopen(OFile,write_mode);
if(isempty(Table)==1)
   fclose all;
   return;
end
peptides=length(Table(:,1));
column=length(Table(1,:));
if(isempty(header)==0)
    fprintf(file_indx,[header,'\n']);  % write names of coloumns
end
for j=1:1:peptides   %for all best signal peptides
    for c=1:1:column-1
        fprintf(file_indx,'%s\t',num2str(Table{j,c}));  %% file save best SP in an index
    end
    fprintf(file_indx,'%s\n',num2str(Table{j,column}));  %% file save best SP in an index
end
fclose all;
end