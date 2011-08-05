function []=CompareExpTheoPeptides(TheoIF,ExpIF,NotDetIF,col,header,sel)

%% Compare theoretical Peptides with Experimental
%%  INPUT
%%      DatabaseIF: File path
%%      col:        [X Y Z K]
%%                  X,Y: column where ID is in Theoretical and Experimental File
%%                  Z,K: column where aa sequence is in Theoretical and Experimental File
%%      header:     [A B] 1 if the corresponding file contains header file contains header
%%                  0 otherwise
%%  OUTPUT
%%      Table:      a table containing all ID (one per line) and the number
%%                  of stances per ID

% ID | Sequence |LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity|Serial
[namefile dir]=IsolateFileName({NotDetIF});
[namefile_th dir]=IsolateFileName({TheoIF});
[namefile_exp dir]=IsolateFileName({ExpIF});
PepExp_if=['Data/TrypticPeptides/PEPexp_',namefile_exp{1},'.txt'];
PepTheo_if=['Data/TrypticPeptides/PEPtheo_',namefile_th{1},'.txt'];
PEPNonDet_if=['Data/TrypticPeptides/PEPNonDet_',namefile{1},'.txt'];
CountExp_if=['Data/TrypticPeptides/COUNTexp_',namefile_exp{1},'.txt'];
CountTheo_if=['Data/TrypticPeptides/COUNTtheo_',namefile_th{1},'.txt'];
CountNonDet_if=['Data/TrypticPeptides/COUNTNonDet_',namefile{1},'.txt'];
    
f=[fopen(PepTheo_if) fopen(PepExp_if) fopen(PEPNonDet_if)];sel=sel.*(f>0);
fclose all;
if (sel(1)==1)
    [A B C D E F G H I] = textread(PepTheo_if,'%s %s %f %f %f %f %f %f %f',-1,'delimiter','\t');
    UNIQUE_PEP_Theo=[A Double2CellTable(B) Double2CellTable(C) Double2CellTable(D) Double2CellTable(E) Double2CellTable(F) Double2CellTable(G) Double2CellTable(H) Double2CellTable(I)];
    [UniqueIDs_Theo(:,1) UniqueIDs_Theo(:,2)] = textread(CountTheo_if,'%s %s',-1,'delimiter','\t');
else
    [UniqueIDs_Theo UNIQUE_PEP_Theo]=CountLineSameID(TheoIF,[col(3) col(4)],header(1));
    FileWriteTable(PepTheo_if,UNIQUE_PEP_Theo,[],'w');
    FileWriteTable(CountTheo_if,UniqueIDs_Theo,[],'w');    
end
if (sel(2)==1)
    [A B C D E F G H I]= textread(PepExp_if,'%s %s %f %f %f %f %f %f %f',-1,'delimiter','\t');
    UNIQUE_PEP_Exp=[A Double2CellTable(B) Double2CellTable(C) Double2CellTable(D) Double2CellTable(E) Double2CellTable(F) Double2CellTable(G) Double2CellTable(H) Double2CellTable(I)];
    [UniqueIDs_Exp(:,1) UniqueIDs_Exp(:,2)]= textread(CountExp_if,'%s %s',-1,'delimiter','\t');
else
    [UniqueIDs_Exp UNIQUE_PEP_Exp]=CountLineSameID(ExpIF,[col(1) col(2)],header(2));
    FileWriteTable(PepExp_if,UNIQUE_PEP_Exp,[],'w');
    FileWriteTable(CountExp_if,UniqueIDs_Exp,[],'w');    
end
if (sel(3)==1)
    [A B C D E F G H I]= textread(PEPNonDet_if,'%s %s %f %f %f %f %f %f %f',-1,'delimiter','\t');
    UNIQUE_PEP_NonDet=[A Double2CellTable(B) Double2CellTable(C) Double2CellTable(D) Double2CellTable(E) Double2CellTable(F) Double2CellTable(G) Double2CellTable(H) Double2CellTable(I)];
    [UniqueIDs_NonDet(:,1) UniqueIDs_NonDet(:,2)]= textread(CountNonDet_if,'%s %s',-1,'delimiter','\t');            
else
    [UniqueIDs_NonDet UNIQUE_PEP_NonDet]=CountLineSameID(NotDetIF,[col(5) col(6)],header(3));
    FileWriteTable(PEPNonDet_if,UNIQUE_PEP_NonDet,[],'w');
    FileWriteTable(CountNonDet_if,UniqueIDs_NonDet,[],'w');    
end

% figure(2);scatter(MOLWtheo(:,1),MOLWtheo(:,2),'k','Marker','+');hold on;
% figure(2);scatter(MOLWexp(:,1),MOLWexp(:,2),'r','Marker','+');hold off;
% xlabel('Molecular Weight (kD)');
% ylabel('Length (aas)');
% figure(3);scatter(MOLWtheo(:,1),MOLWtheo(:,3),'k','Marker','+');hold on;
% figure(3);scatter(MOLWexp(:,1),MOLWexp(:,3),'r','Marker','+');hold off;
% xlabel('Molecular Weight (kD)');
% ylabel('GRAVY Index');

% ID | Sequence |LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity

% peptides=size(UNIQUE_PEP_Exp,1);
% index_of_notdet=zeros(size(UNIQUE_PEP_Theo,1),1);
% index_of_notdet=logical(index_of_notdet);
% TheoSeq=CellTable2StrTable(UNIQUE_PEP_Theo(:,2));
% for i=1:peptides
%     found=strcmpi(TheoSeq,UNIQUE_PEP_Exp{i,2});   %% compare all IDs with cur ID to see which are detected
%     index_of_notdet=(index_of_notdet | (not(index_of_notdet) & found));
% end

[Not_FOUND_m2z Not_FOUND_byfeatures One_Pep_One_Charge]=Plot3D(UNIQUE_PEP_Exp,UNIQUE_PEP_Theo,UNIQUE_PEP_NonDet,25);  %max(CellTable2Double(UniqueIDs_Theo(:,2)))
% [Not_FOUND_m2z Not_FOUND_byfeatures One_Pep_One_Charge]=Plot3D(UNIQUE_PEP_Exp,UNIQUE_PEP_Theo(not(index_of_notdet),:),UNIQUE_PEP_NonDet);

proteins=size(UniqueIDs_Exp,1);
SummarizeTable=[UniqueIDs_Theo cell(size(UniqueIDs_Theo,1),1)];

for i=1:proteins
    found=find(strcmpi(UniqueIDs_Theo,UniqueIDs_Exp(i,1)),1);   %% compare all IDs with cur ID
    if(isempty(found)==0)
        SummarizeTable(found,3)=UniqueIDs_Exp(i,2);
    end
end
THEOR=CellTable2Double(SummarizeTable(:,2));
EXP=CellTable2Double(SummarizeTable(:,3));
DSummarizeTable=( EXP./THEOR ) *4;
DSummarizeTable(DSummarizeTable>4)=4;
peptidenum=1:5:4;
TEMP=Freqcalc(DSummarizeTable,peptidenum);
[h1]=FigureLegends(peptidenum,TEMP',1,'Tryptic Peptides Detected(Percent over theoretical)','Percent of Proteins','Theoretical VS Experimental Tryptic Peptides',{'Theoretical' 'Experimental'},'b',{'-','';':','o'});
folder=regexp(ExpIF,'[/]');
dot=regexp(ExpIF,'[.]');
if(isempty(folder)==0)
    FileWriteTable([ExpIF(1:folder(end)),ExpIF(folder(end)+1:dot(1)-1),'_COUNT.txt'],SummarizeTable,'ID\tThoretical\tExperimental','w');
else
    FileWriteTable([ExpIF(1:dot(1)-1),'_COUNT.txt'],SummarizeTable,'ID\tThoretical\tExperimental','w');
end


end