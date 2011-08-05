function []=FastaSubSelection(INames,IFasta,idtype,col)
%% INPUT:
%% 	INames      file with a list of accessions or mnemonics or gene names
%% 	IFasta      fasta file in which you search for the above names
%% 	idtype      what type of id is in file "INames"
%% 					'mn':		for mnemonics
%% 					'acc':	for accession IDs
%% 					'gn':		for gene names
%%  col         [X Y] X coloumn were ID is Y if header exists in <INames>
%%              file
%% 	wr_acc      1 for printing accession id in separate colouns in output file
%%
%% OUTPUT:
%%  <[IFasta.INames].fasta>     Fasta file with IDs from "INames" file
%%  <[IFasta.INames]_.fasta>    Fasta file with IDs that are not in "INames" but in fasta file (negative set)
%%  <[IFasta]_mn.txt>           file with found IDs
%%  <NotFoundBySequence.txt>     a file that contains IDs not found in specified
%%                              fasta "IFasta"

found=regexp(IFasta,'[/.]');
if(isempty(found)==0)
    dirfasta=IFasta(1:found(end));
else
    dirfasta=[];
end
found=regexp(INames,'[/]');name=regexp(INames,'[.]');
if(isempty(found)==0)
    dirdataset=INames(1:found(end));
    fileName=INames(found(end)+1:name(end)-1);
else
    dirdataset=[];
    fileName=INames(1:name(end)-1);
end

[HEAD SEQUENCE]=FastaRead(IFasta);
% [GnNames, TotLen, SEQ] = textread(INames,'%s %d %s',-1,'delimiter','\t');SPLen=zeros(length(SEQ));  %%% read dataset file
% [GnNames, TotLen] = textread(INames,'%s %d',-1,'delimiter','\t');SEQ=cell(1,length(GnNames));SPLen=zeros(length(SEQ));  %%% read dataset file
% [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(INames,'%s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file
% [GnNames, Descr,TotLen, MW, SPLen, SEQ] = textread(INames,'%s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file
% [GnNames, Accessions, Descr, Comments, TotLen, MW, SPLen, SEQ] = textread(INames,'%s %s %s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file
% [GnNames, Comments] = textread(INames,'%s %s',-1,'delimiter','\t');
%% Read MUTANT
[Names] = ReadTable(INames);
GnNames=Names(col(2)+1:end,col(1));
Discr=cell(1,length(GnNames));SEQ=cell(1,length(GnNames));SPLen=zeros(1,length(SEQ));  %%% read dataset file
%% REGEXP for Mnemonic, Gene Name, Accession Number
mnem_expr='[A-Za-z0-9]{3,6}_ECOL[A-Z]';
% mnem_expr='([A-Za-z0-9]{3,6}_MYCTU) | ([A-Za-z0-9]{3,6}_ECOL[A-Z])';
acc_expr='[A-Z][0-9A-Za-z]{5}(-[0-9])*[|]';
gn_exps='=[a-z]{3}[A-Z0-9]*';
[GNs_st GNs_en GNs_ext GNs_]=regexp(HEAD,gn_exps); % Gene Names
[GNs]=CellTable2StrTable_(GNs_,'gn');
[Acc_st Acc_en Acc_ext Acc_]=regexp(HEAD,acc_expr); % Gene Accessions
[Acc]=CellTable2StrTable_(Acc_,'acc');
[MNs_st MNs_en MNs_ext MNs_]=regexp(HEAD,mnem_expr); % Gene Mnemonics
[MNs]=CellTable2StrTable_(MNs_,'mn');
switch idtype
    case 'mn'
        Comp=MNs;expr=mnem_expr;
    case 'gn'
        Comp=GNs;expr=gn_exps;
    case 'acc'
        Comp=Acc;expr=acc_expr;
    otherwise
        display('=>Wrong type ID mn/acc <FastaSubSelection.m>');
        return;
end
%% Init File with not found proteins
WriteMnemonics([dirdataset,'NotFoundBySequence.txt'],{'hi'},'w');
%% Init loop
allheaders=length(HEAD);
remaining=ones(1,allheaders);
headers=length(GnNames);
HEAD_NEW=cell(1,headers);
SEQ_NEW=cell(1,headers);
MNEM=cell(1,headers);
ACC=cell(1,headers);
TotLen=zeros(1,headers);
MW=zeros(1,headers);
len=11;indx_seq_sub=[];
sub_SEQUENCE=SubStrings(SEQUENCE,len);
for i=1:1:headers
    %% Match Sequence or Identifier in fasta
    indx_seq=find(strcmpi(SEQUENCE,SEQ{i}));
    indx_name=find(strcmpi(Comp,GnNames{i}));indx=indx_name;
    if(isempty(SEQ{i})==0)
        indx_seq_sub=find(strcmpi(sub_SEQUENCE,SEQ{i}(1:len)));
    end
    if(isempty(indx_seq)==0)
        [st en ext match]=regexp(HEAD{indx_seq},expr);indx=indx_seq;
    elseif(isempty(indx_name)==0)
        if(length(indx_name)>1)
            display([GnNames{i},' found twice']);
        end
        [st en ext match]=regexp(HEAD{indx_name(1)},expr);indx=indx_name(1);
    elseif(isempty(indx_seq_sub)==0)
        [st en ext match]=regexp(HEAD{indx_seq_sub},expr);indx=indx_seq_sub;
    end
    
    %% In case there are amino acids 'X,U' replace with 'A'
    if(isempty(indx)==0)
        remaining(indx)=0;
        [st]=regexp(SEQUENCE{indx},'[UX]+');
        if(isempty(st)==0)
            SEQUENCE{indx}(st)='A';
        end
        MNEM{i}=MNs{indx};          % Mnemonic of found
        ACC{i}=Acc{indx};           % Accession of found
        HEAD_NEW{i}=HEAD{indx};     % header of fasta
        SEQ_NEW{i}=SEQUENCE{indx};  % Sequence of fasta
        TotLen(i)=length(SEQUENCE{indx});       % protein length
        MW(i)=1;%molweight(SEQUENCE{indx});
    else
        display(['Not found: ',GnNames{i}]);
        WriteMnemonics([dirdataset,'NotFoundBySequence.txt'],[{GnNames{i}}],'a');
    end
end

empty_indx=ones(1,headers);
for i=1:headers
    if(isempty(MNEM{i}))
        empty_indx(i)=0;
    end
end
empty_indx=logical(empty_indx);

WriteMnemonics([dirdataset,fileName,'_Mnemonic.txt'],{MNs{logical(remaining)}},'w');
WTable=[MNEM' ACC' HEAD_NEW' Double2CellTable(TotLen') Double2CellTable(round(MW)') Double2CellTable(SPLen') SEQ_NEW'];
FileWriteTable([dirdataset,fileName,'_mn.txt'],WTable(empty_indx,:),[],'w');
FastaWrite([dirfasta,fileName,'.fasta'],HEAD_NEW,SEQ_NEW,70,'w');
FastaWrite([dirfasta,fileName,'_.fasta'],{HEAD{logical(remaining)}},{SEQUENCE{logical(remaining)}},70,'w');

end

function []=WriteMnemonics(oFile,MNEM,open)
s=fopen(oFile,open);  %% file save group of proteins with specific best SP
for i=1:1:length(MNEM)
    fprintf(s,[MNEM{i},'\n']);
end
fclose(s);
end

function [SubStrings]=SubStrings(StringList,len)
STRINGS=length(StringList);
SubStrings=cell(STRINGS,1);
for i=1:1:STRINGS
    SubStrings{i}=StringList{i}(1:len);
end

end

%% Transform regexp output cell x cell table to cell of stings
function [Table]=CellTable2StrTable_(CellTable,idtype)
tokens=length(CellTable);
Table=cell(tokens,1);

for i=1:1:tokens
   if(isempty(CellTable{i}))
       Table(i)={'0'};
   else
       switch idtype
           case 'mn'
               Table{i}=CellTable{i}{1,1};
           case 'gn'
               Table{i}=CellTable{i}{1,1}(2:end);
           case 'acc'
               Table{i}=CellTable{i}{1,1}(1:end-1);
           otherwise
               display('=>Wrong type ID mn/acc <FastaSubSelection.m>');
               return;
       end
   end
end
end
