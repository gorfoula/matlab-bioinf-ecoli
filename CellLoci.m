function []=CellLoci(ifCatgCorresp,ifProteomeDB,ifIdentifiedList,comment_merge,catg,indexing,header)
%% INPUT
%%  ifCatgCorresp:      File which contains table => [ID|Category Name]
%%  ifProteomeDB:       Data file (column of MW is now Catg ID)
%%  ifIdentifiedList:   List of names to be assigned to a category
%%  comment_merge:      1 to merge comment from 2 files
%%  indexing:           [A B] <int> 
%%                      where A is coloumn number where ID (e.g Accession)
%%                      is in <ifProteomeDB> file
%%                      and B the coloumn number where the corresponding ID
%%                      is in <ifIdentifiedList> file
%%  catg:               Coloumn number in <ifProteomeDB> file were Category
%%                      annotation is
%%  header:             [X y] header existance in <ifProteomeDB> and
%%                      <ifIdentifiedList>
%%                      1   if header exists
%%                      0   no header    

%% READ Files
OFile=[ifIdentifiedList,'_OUT.txt'];
[IDs, Catg] = textread(ifCatgCorresp,'%s %s',-1,'delimiter','\t');  %%% read dataset file

ProteomeDB_=ReadTable(ifProteomeDB);
IdentifiedList_=ReadTable(ifIdentifiedList);

ProteomeDB=ProteomeDB_(header(1)+1:end,:);
IdentifiedList=IdentifiedList_(header(2)+1:end,:);

%% Match Mnemonics in List

[pos_db not_found doubles found]=MatchList(IdentifiedList(:,indexing(2)),ProteomeDB(:,indexing(1)));
pos_db=pos_db(pos_db>0); % were in database they were found
not_found=not_found(not_found>0);   % those not found in the database
found=found(found>0);  % index on the Scanned List
doubles=doubles(doubles>0);  % doubles index on database

if(isempty(doubles)==0)
    FileWriteTable([OFile,'_doubles.txt'],ProteomeDB(doubles,:),header,'w');
end
FileWriteTable(OFile,ProteomeDB(pos_db,:),header,'w');
if(isempty(comment_merge)==0)
    Comb=cellfun(@(a,b) [a,b],ProteomeDB(pos_db,comment_merge(1)),IdentifiedList(found,comment_merge(2)),'uni',false);
    FileWriteTable([OFile,'_updated.txt'],[ProteomeDB(pos_db,1:2) Comb ProteomeDB(pos_db,4:end)],header_new,'w');
else
    FileWriteTable([OFile,'_updated.txt'],ProteomeDB(pos_db,:),header,'w');
end
if(isempty(not_found)==0)
    FileWriteTable([OFile,'_updated.txt'],ProteomeDB(not_found,:),header,'a');
end
DrawPie(ProteomeDB(pos_db,catg),ProteomeDB(:,catg),IDs,Catg);

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
    if(i==2623)
        display('hey');
    end
    indx_found=find(strcmpi(GnNames,Mnemonic_{i}));
    if(isempty(indx_found)==0)          % if proteins found in total list
        if(sum(pos_found==indx_found(end))==0)  % if position found is not already in index list, add it
            pos_found(i)=indx_found(end);
            not_found(pos_found(i))=0;
            found(i)=i;
        else
            display(['Double in scanned List: ',Mnemonic_{i}])
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
proteins=length(DataBase);
dataset_proteins=length(CATGs);

for i=1:1:categories-1
    indx_found=find(strcmpi(CATGs,IDs(i)));
    indx_database=find(strcmpi(DataBase,IDs(i))); %% found in whole database
%     Nfound(i)=(length(indx_found)./proteins)*100;   %% Number fond from each category
    Nfound(i)=length(indx_found);   %% Number fond from each category
%     Catg(i)={[Catg{i},' (',num2str(Nfound(i)),')']};
    Catg(i)={[Catg{i},' (',num2str(ceil(Nfound(i)*100./dataset_proteins)),'%)']};
    PercentOverAll(i)=Nfound(i)./length(indx_database);
end

explode = zeros(categories,1);
explode(1)=1;
figure(1);subplot(2,2,3);
h=pie3(Nfound(Nfound>0),explode(Nfound>0),Catg(Nfound>0));ColorPie(h,Nfound(Nfound>0));
figure(1);subplot(2,2,[1 2]);
% bar(1:1:categories,PercentOverAll,'y');hold on;
h = stem(1:1:categories,PercentOverAll,'LineWidth',5,'Color','k');
set(h,'Marker','none');
set(gca,'XTick',1:1:categories);
set(gca,'XTickLabel',IDs);

for i=1:1:categories
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
Colors=RainbowColor(100);
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

