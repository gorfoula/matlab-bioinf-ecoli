function [TABLE index files print_cols]=ExludeCyto(Files,file_ref,rerun)
%% 
%%  INPUT:
%%      Files       File which contains the file list (Input of <JoinMSMSResults>)
%%                  Database file must be fist in the file list
%%                  Database file is used to exclude older version
%%                  Accessions
%%      file_ref    Which file is our file in the file list (e.g 2 if it is the second in order)
%%                  this is used to calculate novel proteins
%%      rerun       Rerun merging and don't load previous files
%%      e.g. 
%%      JoinOtherStudiesScript('Notepad Lists/OthersFiles.txt',7);


%%% LOAD FILE LIST   )))))))))))))))))))))))
[FileNames] = ReadTable(Files,'\n');
files=size(FileNames,1)-1;
printarray=CellTable2Double(FileNames(:,4:end));
print_cols=[sum(printarray>0,2) printarray];
%%% DISCARD FILE PATHWAY   )))))))))))))))))))))))
[filename directory]=IsolateFileName({Files});filename=filename{1,1};directory=directory{1,1};
%%% DEFINE OUTPUT FILES   ))))))))))))))))))
Out_File=[directory,'(',filename,')','_COMPARISON.txt'];
Out_File_CellEnv=[directory,'(',filename,')','_COMPARISON_CELL_ENV.txt'];
Out_File_Cyto=[directory,'(',filename,')','_COMPARISON_CYTO.txt'];
Out_File_CellEnv_mark=[directory,'(',filename,')','_COMPARISON_CELL_ENV_MARKED.txt'];
%%% LOAD PREVIOUS FILES OR RERUN MERGE   ))))))
file_index=fopen(Out_File_CellEnv);
file_index_cyto=fopen(Out_File_Cyto);
if(file_index>0 && file_index_cyto>0 && rerun==0)
    fclose(file_index);
    fclose(file_index_cyto);
    [CELLENV] = ReadTable(Out_File_CellEnv,'\n');
    CELLENV=[CELLENV(1,:);CELLENV(3:end,:)];
    [CYTO] = ReadTable(Out_File_Cyto,'\n');
    TABLE=[CELLENV;CYTO];
    index=[ones(size(CELLENV,1)-1,1);zeros(size(CYTO,1),1)];
    display(['CAUTION: File <',Out_File_CellEnv,'> is loaded!']);
    return;
end
%%% MERGE FILES  )))))))))))))))))))))))
[list,Header,files,discr]=JoinMSMSResults(Files,1);
Check_Table=CellTable2Double(list(:,2:end));
%% EXCLUDE ACCESSIONS NOT IN CURENT DB ))))))))))))))))))))
In_curated_DB=(Check_Table(:,1)==1);% & (JustNumbers_comp(:,files+1)>1);
list=list(In_curated_DB,:);
discr=discr(In_curated_DB,:);
Check_Table=Check_Table(In_curated_DB,:);
%%% REMOVE COLOUMNS OF DB  )))))))))))))))))))
Check_Table(:,files+1)=Check_Table(:,files+1)-1;
Check_Table=[Check_Table(:,2:files+1) Check_Table(:,files+3:end)];
Header=[Header(1);Header(3:files+3);Header(2*files+3+1:end)];
files=files-1;
%%%  HOW MANY PROTEINS FROM EACH FILE   ))))))))))))))))))
each_col_sum=zeros(1,files+1);
each_col_sum(1:files)=sum(Check_Table(:,1:files));
each_col_sum(files+1)=sum(Check_Table(:,files+1)==1 & Check_Table(:,file_ref)==1);
Count_str=Double2CellTable(each_col_sum);
Count_str=[{'#Proteins'} Count_str];
Count_str(end)={['Novel: ',num2str(Count_str{end})]};
%%%  VENN )))))))))))))
if(files<=3 && files>1)
    VenDiagram(Check_Table(:,1:files),Header(2:files+1));
    title('K12 TOTAL PROTEOME');
end
%%% IN HOW MANY STUDIES EACH PROTEIN IS FOUND  )))))))))))))
figure;[freq x]=hist(Check_Table(:,files+1),0:1:max(Check_Table(:,files+1)));bar(x,(freq*100)./size(Check_Table,1));
title('K12 TOTAL PROTEOME');
xlabel('# of studies');
ylabel('# of proteins');
%%%%  PRINT ALL PROTEINS  )))))))))))))))))
New_list=[list(:,1) Double2CellTable(Check_Table(:,1:files+1)) discr];
FileWriteTable(Out_File,Header',[],'w');
FileWriteTable(Out_File,Count_str,[],'a');
FileWriteTable(Out_File,New_list,[],'a');
%%% WHICH ARE CEP?  )))))))))))))))))))
index_CellEnv=not( strcmp(New_list(:,files+3),'A') | strcmp(New_list(:,files+3),'A_trl') );
%%%  VENN )))))))))))))
if(files<=3 && files>1)
    VenDiagram(Check_Table(index_CellEnv,1:files),Header(2:files+1));
    title('K12 CELL ENVELOPE');
end
%%%  COUNT HOW MANY PROTEINS FOUND IN EACH FILE )))))))))))))
[Count_str]=CountOccurances(Check_Table,files,file_ref);
%%%  HOW MANY OF THE CEP? )))))))))))))
[Count_str_cellenv]=CountOccurances(Check_Table(index_CellEnv,:),files,file_ref);
%%%  PRINT CEP PROTEINS   )))))))))))))
FileWriteTable(Out_File_CellEnv,Header',[],'w');
FileWriteTable(Out_File_CellEnv,Count_str_cellenv,[],'a');
FileWriteTable(Out_File_CellEnv,New_list(index_CellEnv,:),[],'a');
%%%  PRINT ALL PROTEINS FIRST CEP PROTEINS   )))))))))))))
CELLENV=[Header';New_list(index_CellEnv,:)];
CYTO=[Header';New_list(not(index_CellEnv),:)];
FileWriteTable(Out_File_CellEnv_mark,Header',[],'w');
FileWriteTable(Out_File_CellEnv_mark,Count_str,[],'a');
FileWriteTable(Out_File_CellEnv_mark,Count_str_cellenv,[],'a');
FileWriteTable(Out_File_CellEnv_mark,CELLENV,[],'a');
FileWriteTable(Out_File_CellEnv_mark,CYTO,[],'a');
%%%  PRINT CYTOPLASMIC PROTEINS  )))))))))))))
FileWriteTable(Out_File_Cyto,New_list(not(index_CellEnv),:),[],'w');
%%% IN HOW MANY STUDIES EACH CEP PROTEIN IS FOUND  )))))))))))))
figure;[freq x]=hist(Check_Table(index_CellEnv,files+1),0:1:files);bar(x,(freq*100)./sum(index_CellEnv));
xlabel('# of studies');
ylabel('# of proteins');
title('K12 CELL ENVELOPE');
display('End of Exclude Cyto!');

TABLE=[CELLENV;CYTO];
index=[ones(size(CELLENV,2),1);zeros(size(CYTO,2),1)];
end

function [Count_str]=CountOccurances(MergeTable,files,file_ref)

each_col_sum=zeros(1,files+1);
each_col_sum(1:files)=sum(MergeTable(:,1:files));
each_col_sum(files+1)=sum(MergeTable(:,files+1)==1 & MergeTable(:,file_ref)==1);
Count_str=Double2CellTable(each_col_sum);
Count_str=[{'# Proteins =>'} Count_str];
Count_str(end)={['Novel: ',num2str(Count_str{end})]};

end
