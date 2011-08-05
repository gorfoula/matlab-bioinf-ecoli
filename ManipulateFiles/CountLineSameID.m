function [UniqueIDs UNIQUE_PEP]=CountLineSameID(DatabaseIF,col,header)

%% Counts instances (lines) per ID
%%  INPUT
%%      DatabaseIF: File path
%%      col:        [A B]   A: which coloumn contains the ID
%%                          B: which coloumn contains the Sequence
%%      header:     1 if the file contains header
%%                  0 otherwise
%%  OUTPUT
%%      Table:      a table containing all ID (one per line) and the number
%%                  of instances per ID
FileReadTable=ReadTable(DatabaseIF);
lines=size(FileReadTable,1);
IDs=FileReadTable(header+1:end,col(1)); %% ID lines
PeptideSequence=FileReadTable(header+1:end,col(2)); %% Peptide Sequence lines

UniqueIDs=[];
%%%%% OUTPUT FILE  %%%%%
folder=regexp(DatabaseIF,'[/]');
dot=regexp(DatabaseIF,'[.]');
if(isempty(folder)==0)
    FileName=[DatabaseIF(1:folder(end)),DatabaseIF(folder(end)+1:dot(1)-1),'_peptides.txt'];
else
    FileName=[DatabaseIF(1:dot(1)-1),'_peptides.txt'];
end
FileWriteTable(FileName,{'BEGIN:'},[],'w');
UNIQUE_PEP=[];SerialNum=[];
display(['File: <',DatabaseIF,'>']);
count=1;
while (isempty(IDs)==0)
    cur_ID=IDs(1);
    matches=strcmpi(IDs,cur_ID);   %% compare all IDs with cur ID
    [AllCleaved]=TotalCleave(PeptideSequence(matches));
    [counts UniqueCleaved]=UniquePeptides(AllCleaved,cur_ID{1},FileName);
    IDs=IDs(not(matches));       %% exlude current ID
    PeptideSequence=PeptideSequence(not(matches));
    UniqueIDs=[UniqueIDs;[cur_ID num2str(counts)]];
    SerialNum=[SerialNum;ones(size(UniqueCleaved,1),1)*count];
    UNIQUE_PEP=[UNIQUE_PEP;UniqueCleaved];
    count=count+1;
end
UNIQUE_PEP=[UNIQUE_PEP Double2CellTable(SerialNum)];

end

function [COUNT UniqueCleaved]=UniquePeptides(Sequence,ID,FileName)
COUNT=0;
FileWriteTable(FileName,{ID},[],'a');
UniqueCleaved=[];
while (isempty(Sequence)==0)
    cur_Seq=Sequence{1};
    [copies EXACT ALL]=StrFindArrays(Sequence,cur_Seq);
    cur_Seq(strfind(cur_Seq,'X'))='A';
    cur_Seq(strfind(cur_Seq,'U'))='A';
    [PROPERTIES]=MolWeightCell(cur_Seq);  % LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity
%     if((PROPERTIES{2}>800) & (PROPERTIES{3}>=2) & (PROPERTIES{4}<1800))
        COUNT=COUNT+1;
        FileWriteTable(FileName,{cur_Seq},[],'a');
        UniqueCleaved=[UniqueCleaved;[ID Sequence(1) PROPERTIES]]; % ID | Sequence |LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity
%     else
%         FileWriteTable([FileName,'_outInstrRange.txt'],[ID Sequence(1) PROPERTIES],[],'a');
%     end
    Sequence=Sequence(not(EXACT));
end
    
end

function [AllCleaved]=TotalCleave(Sequence)
AllCleaved=Sequence;
if(isempty(Sequence))
    return;
else
    count=1;
    for i=1:size(Sequence,1)
%         cur_Seq=AllCleaved{i}
        [start_idx, end_idx, extents, matches]=regexp(AllCleaved{count},'[QWETYIASDFGHLPCVNM]*([RK]{0,1})([P][QWETYIASDFGHPLCVNM]*[RK]{0,1})*');
        No_peptides=size(matches,2);
        if(No_peptides>1)
            AllCleaved=[AllCleaved(1:count-1);matches';AllCleaved(count+1:end)];count=count+No_peptides;
        else
            count=count+1;
        end
    end
end
end


function [PROPERTIES]=MolWeightCell(Sequence) % LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity

Sequence(strfind(Sequence,'X'))='A';
Sequence(strfind(Sequence,'U'))='A';
PEP_LEN=size(Sequence,2);
[Indx SelSeq]=AARepresentation('int',1,PEP_LEN,{Sequence});
[m2z PEP_MASS PEP_CHARGE]=Mass2Charge(Sequence);
[Spear PEP_PHO]=Hydrophobicity(Indx,'K');
[PolMatrix PEP_POL]=Polarity(Indx,'Z');
PROPERTIES=[{PEP_LEN} {PEP_MASS} {PEP_CHARGE} {m2z} {PEP_PHO./PEP_LEN} {PEP_POL./PEP_LEN}];

end