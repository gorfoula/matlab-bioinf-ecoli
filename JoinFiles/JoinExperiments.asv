function []=JoinExperiments(file,index)

%%  INPUT:
%%      file:   input file of 
%%      index:  [A B]
%%              A:  how many lines are occupied by the header (0 if no header)
%%              B:  coloumn with experiment condition 

[INTABLE] = ReadTable(file,'\n');
[namefile dir]=IsolateFileName({file});
outfile=[dir{1},'Merge_',namefile{1},'.txt'];
INTABLE_noheader=INTABLE(index(1)+1:end,:);
header=INTABLE(index(1),2:end);

[TABLE Address]=Douplicates(file,index,INTABLE);
FileNames=cell(max(Address),3);

for i=unique(Address)'
    cur_INTABLE=INTABLE_noheader(Address==i,2:end);
    FileNames{i,1}=[dir{1},TABLE{i,1},'.txt']; %% file path
    FileNames(i,2:3)=[{'1'} {'1'}];
    FileWriteTable(FileNames{i,1},[header;cur_INTABLE],[],'w');
end

FileWriteTable(outfile,FileNames,[],'w');

[list,Header,files,discr]=JoinMSMSResults(outfile,1);%% out file

list_=CellTable2Double(list(:,2:files+1));
sum_=sum(list_,2);
uique_=zeros(files,1);
for f=1:files
    uique_(f)=sum((list_(:,f).*sum_)==1);
end
uique_
outfile2=[dir{1},'Unique',namefile{1},'.txt'];
FileWriteTable(outfile2,[Header(1:files+2)';[{''} Double2CellTable(uique_)' {''}];list(:,1:files+2)],[],'w');
end