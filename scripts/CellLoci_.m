function []=CellLoci(IFCatgID,DatabaseIF,IFScanedList,toolspred,comment_merge)
%% INPUT
%%  IFCatgID:   File which contains table => ID|Category Name
%%  DatabaseIF:     Data file (column of MW is now Catg ID)
%%  IFScanedList:   List of names to be assigned to a category
%%  comment_merge:  1 to merge comment from 2 files

%% READ Files
OFile=[IFScanedList,'_OUT.txt'];
[IDs, Catg] = textread(IFCatgID,'%s %s',-1,'delimiter','\t');  %%% read dataset file
    
    [text] = textread(DatabaseIF,'%s',-1,'delimiter','\n');
    header=text{1};
    [start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text(1:end),'[\t]');
    Table=CellTable2StrTable(Table_);

    [text] = textread(IFScanedList,'%s',-1,'delimiter','\n');
    header_new=text{1};
    [start_idx, end_idx, extents, matches, tokens, names, TableScanned_] = regexp(text(1:end),'[\t]');
    TableScanned=CellTable2StrTable(TableScanned_);

%% Match Mnemonics in List

[pos_db not_found doubles found]=MatchList(TableScanned(:,1),Table(:,2));
pos_db=pos_db(pos_db>0); % were in database they were found
not_found=not_found(not_found>0);   % those not found in the database
found=found(found>0);  % index on the Scanned List
doubles=doubles(doubles>0);  % doubles index on database

if(isempty(doubles)==0)
    FileWriteTable([OFile,'_doubles.txt'],Table(doubles,:),header,'w');
end
FileWriteTable(OFile,Table(pos_db,:),header,'w');
if(comment_merge)
    Comb=cellfun(@(a,b) [a,b],Table(pos_db,3),TableScanned(found,3),'uni',false);
    if(toolspred)
        FileWriteTable([OFile,'_updated.txt'],[Table(pos_db,1:6) TableScanned(found,7:end)],header_new,'w');
    else
        FileWriteTable([OFile,'_updated.txt'],[Table(pos_db,1:2) Comb Table(pos_db,4:6) TableScanned(found,7:end)],header_new,'w');
    end
else
    FileWriteTable([OFile,'_updated.txt'],Table(pos_db,:),header,'w');
end
if(isempty(not_found)==0)
    FileWriteTable([OFile,'_updated.txt'],Table(not_found,:),header,'a');
end
DrawPie(Table(pos_db,6),Table(:,6),IDs,Catg);

end

%%  Match found proteins with categories
function [pos_found not_found doubles_index found]=MatchList(Mnemonic_,GnNames)
%% INPUT
%%  Mnemonic_:  List of Menmonic to match
%%  GnNames:    All mnemonics
%%  ID:         All mnemonics IDs
%% OUTPUT
%%  pos_found: index of positions found in total file
%%  doubles_index:  double found 

proteinsF=length(Mnemonic_);
not_found=1:1:length(GnNames);
pos_found=zeros(1,proteinsF);
found=zeros(1,proteinsF);
doubles_index=zeros(1,proteinsF);
count_doubles=0;
for i=1:1:proteinsF
    indx_found=find(strcmpi(GnNames,Mnemonic_{i}));
    if(isempty(indx_found)==0)          % if proteins found in total list
        if(sum(pos_found==indx_found(end))==0)  % if position found is not already in index list, add it
            pos_found(i)=indx_found(end);
            not_found(pos_found(i))=0;
            found(i)=i;
        end
        if(length(indx_found)>1)        % if you find that protein is multiple in total list
            if(sum(doubles_index==indx_found(end))==0) % if it is not in the double index list, add it
                doubles_index(i)=indx_found(end);
                count_doubles=count_doubles+length(indx_found)-1;
            end
        end
    else
        display(['Not found: ',Mnemonic_{i}]);
    end
end
count_doubles
end

%% Plot a pie of found
function []=DrawPie(CATGs,DataBase,IDs,Catg)
categories=length(IDs);
Nfound=zeros(categories,1);
PercentOverAll=zeros(categories,1);

for i=1:1:categories
    indx_found=find(strcmpi(CATGs,IDs(i)));
    indx_database=find(strcmpi(DataBase,IDs(i))); %% found in whole database
    Nfound(i)=length(indx_found);   %% Number fond from each category
    Catg(i)={[Catg{i},' (',num2str(Nfound(i)),')']};
    PercentOverAll(i)=Nfound(i);
%     PercentOverAll(i)=Nfound(i)./length(indx_database);
end

explode = zeros(categories,1);
explode(1:2)=1;
figure(1);subplot(2,1,1);
h=pie3(Nfound(Nfound>0),explode(Nfound>0),Catg(Nfound>0));ColorPie(h,Nfound);
figure(1);subplot(2,1,2);
bar(1:1:categories,PercentOverAll,'y');hold on;
for i=1:1:categories
    text(i,0,IDs{i},'HorizontalAlignment','center','VerticalAlignment','bottom');
    if(Nfound(i)>0)
        text(i,PercentOverAll(i),[num2str(round(PercentOverAll(i)*100)),'%'],'Color','r','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
        text(i,PercentOverAll(i),num2str(Nfound(i)),'Color','k','FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end
hold off;
%% Print
Nfound(Nfound>0)
Catg(Nfound>0)
sum(Nfound)

end

function []=ColorPie(h,Nfound)
categories=length(Nfound);
%% Colors
Colors=RainbowColor(50);
col_index=1:floor(length(Colors)./categories):length(Colors);
Colors=Colors(col_index,:);
index=[1 2 3];
for i=1:1:sum((Nfound>0))
    edgecol=Colors(i,:)-(100/256);
    edgecol(edgecol<0)=0;
    set(h(index(1)),'FaceColor',Colors(i,:),'EdgeColor',edgecol);
    set(h(index(2)),'FaceColor',Colors(i,:),'EdgeColor','none');
    set(h(index(3)),'FaceColor',Colors(i,:),'EdgeColor',edgecol);
    index=index+4;
end


end