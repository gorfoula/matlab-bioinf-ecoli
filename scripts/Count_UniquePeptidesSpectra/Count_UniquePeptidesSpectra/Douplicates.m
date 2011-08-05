function [TABLE Address]=Douplicates(file,index,INTABLE)
%% This function counts how many times each ID occures in the file
%% INPUT
%%      file:   path of the given tab delimited table
%%      index:  [X Y]
%%              X 1 if the tab delimited table has one line header
%%                0 no header
%%                n  if header occupies n lines
%%              Y coloumn where ID is (e.g 1 for first coloumn)
%% OUTPUT
%%      <file_douplicates.txt> outputfile in the same directory where given
%%                             file is that for each ID gives the number of
%%                             occurances (number of lines)

if(exist('INTABLE','var')==0)
   [INTABLE] = ReadTable(file,'\n'); 
end
[namefile dir]=IsolateFileName({file});
outfile=[dir{1},namefile{1},'_douplicates.txt'];
IDS_intact=INTABLE(index(1)+1:end,index(2));
IDS=INTABLE(index(1)+1:end,index(2));

lines=size(IDS,1);
TABLE=[];
Address=zeros(size(IDS_intact,1),1);
count=1;
while (isempty(IDS)==0)
    cur_id=IDS(1);
    [index]=strcmp(IDS,cur_id);
    [index_intact]=strcmp(IDS_intact,cur_id);
    TABLE=[TABLE; [cur_id {num2str(sum(index))}] ];
    IDS=IDS(not(index));
    Address(index_intact)=count;
    count=count+1;
end

FileWriteTable(outfile,TABLE,[],'w');

end