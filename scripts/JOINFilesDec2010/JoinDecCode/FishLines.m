function [FOUNDTABLE count_found count_nfound]=FishLines(TheoIF,IDsIF,col,header)

TheoReadTable_=ReadTable(TheoIF);
SearchForReadTable_=ReadTable(IDsIF);

TheoReadTable=TheoReadTable_(header(1)+1:end,:);
SearchForReadTable=SearchForReadTable_(header(2)+1:end,:);
couples=size(col,1);

[filename_TheoIF]=IsolateFileName(TheoIF);
[filename_IDsIF dir]=IsolateFileName(IDsIF);

OF_All=[dir,'(DB)',filename_TheoIF,'(for)',filename_IDsIF,'(',num2str(col(1)),')L.txt'];display(OF_All);
OF_All_db=[dir,'(DB)',filename_TheoIF,'(for)',filename_IDsIF,'(',num2str(col(1)),')DB.txt'];display(OF_All_db);
OF_NotF=[dir,'(DB)',filename_TheoIF,'(for)',filename_IDsIF,'(',num2str(col(1)),')NF.txt'];display(OF_NotF);
OF_F=[dir,'(DB)',filename_TheoIF,'(for)',filename_IDsIF,'(',num2str(col(1)),')F.txt'];display(OF_F);
OF_DUP=[dir,'(DB)',filename_TheoIF,'(for)',filename_IDsIF,'(',num2str(col(1)),')DUP.txt'];display(OF_DUP);

proteins=size(SearchForReadTable,1);
proteins_db=size(TheoReadTable,1);
f_db=zeros(proteins_db,2);
f_indx=zeros(proteins,2);
links=cell(proteins,2);
AllCell=cell(proteins,1);
AllCell_db=cell(proteins_db,1);

[filename_TheoIF]=IsolateFileName(TheoIF);
[filename_IDsIF]=IsolateFileName(IDsIF);
FileWriteTable(OF_DUP,{'Hey'},[],'w');
for c=1:couples
    UniqueIDs_Theo=TheoReadTable(:,col(c,1));
    IDsReadTable=SearchForReadTable(:,col(c,2));
    remain=UniqueIDs_Theo;
    while (CheckEmpty(remain)>0)
        [token remain]=strtok(remain);
        p=1;
        while (p<=proteins)
            [start_idx, end_idx, extents, matches, tokens, names, splits]=regexp(num2str(IDsReadTable{p}),'[ :]');
            for s=1:length(splits)
                Found= strmatch(splits{s}, token,'exact');
                if(isempty(Found)==0  && strcmp(IDsReadTable{p},'')==0)
                    if(length(Found)>1)
                        f_db(Found)=1;
                        Found=Found(1);
                    end
                    f_indx(p,1)=1;
                    f_indx(p,2)=Found;
                    AllCell_db(Found)=MergeColumns([TheoReadTable(Found,:) SearchForReadTable(p,:)],9);
                    AllCell(p)=MergeColumns([TheoReadTable(Found,:) SearchForReadTable(p,:)],9);
                else
                    if(f_indx(p,1)==0)
%                         links{p}=['=HYPERLINK("http://www.uniprot.org/uniprot/?query="&','A',num2str(p),'&"+ECOLI&sort=score")'];
                        links{p}=['=HYPERLINK("http://www.uniprot.org/uniprot/?query=',IDsReadTable{p},'&sort=score")'];
                        AllCell(p)=MergeColumns([{[9,links{p}]} SearchForReadTable(p,:)],9);
                    end                    
                end
            end
            p=p+1;
        end
    end
end
FileWriteTable(OF_DUP,TheoReadTable(logical(f_db),:),[],'w');
FileWriteTable(OF_F,[TheoReadTable(f_indx(logical(f_indx(:,1)),2),:) SearchForReadTable(logical(f_indx(:,1)),:)],[],'w');
FileWriteTable(OF_NotF,[IDsReadTable(not(f_indx(:,1)),:) SearchForReadTable(not(f_indx(:,1)),:) links(not(f_indx(:,1)))],[],'w');
FileWriteTable(OF_All,AllCell,[],'w');
FileWriteTable(OF_All_db,AllCell_db,[],'w');

FOUNDTABLE=[ [TheoReadTable_(1,:);TheoReadTable(f_indx(logical(f_indx(:,1)),2),:)] [SearchForReadTable_(1,:);SearchForReadTable(logical(f_indx(:,1)),:)] ];
count_found=sum(f_indx(:,1));
count_nfound=sum(f_indx(:,1)==0);
display(['All:    ',num2str(count_nfound+count_found),' Found: ',num2str(count_found),'   NotFound:   ',num2str(count_nfound)]);
display(proteins);
end

function [FirstFound]=CheckEmpty(CellArray)

rows=length(CellArray);
FirstFound=0;
for i=1:rows
    if(isempty(CellArray{i})==0)
        FirstFound=i;
        return;
    end
end

end