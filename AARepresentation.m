function [AAs SelSeq SelSeqend]=AARepresentation(type,start,len,SeqList)
% TYPE      type of representation
%               - 'char':   representation with characters, e.q. 'M' at pos
%               - 'bin':    binary, values [0,1]
%                           for each amino acid a vector of 20 values is allocated, 
%                           1 at the position of the amino according to the representation (below)
%                           Methionine(1) so 1000000000..... for the fist AA
%               -'int':     Integer (1) for Mthionine, (2) Alanine (Representation below)
%LEN        number of Amino acids taken from the total sequence
%SEQLIST    matrix (peptides X 1) with the total sequence of peptides

%--------------------------- Hydrophobic/NonPolar  (7+1)
%Methionine     M   (1)
%Alanine        A   (2)
%Valine         V   (3)
%Leucine        L   (4)
%Isoleucine     I   (5)
%Proline        P   (6)
%Phenylalanine  F   (7)
%Tryptophan     W   (8) --- Hydrophobic/Polar   (7+1)
%--------------------------
%Glycine        G   (9)
%Serine         S   (10)
%Cysteine       C   (11)
%Asparagine     N   (12)
%Glutamine      Q   (13)
%Tyrosine       Y   (14)
%Threonine      T   (15)
%---------------------------Polar/positive/Basic
%Lysine         K   (16)
%Arginine       R   (17)
%Histidine      H   (18)
%---------------------------Polar/negative/acidic
%Aspartic Acid  D   (19)
%Glutamic Acid  E   (20)
AMINOS=20;
[Row,Col,Zax]=size(SeqList);
SelSeq=cell(Row,1);
SelSeqend=cell(Row,1);
%%%%%%%%%%%%%%%%%%%%%%%%   binary Representation   %%%%%%%%%%%%%%%%%%%%%%
if(strcmp(type,'bin'))
    AAs=zeros(Row,AMINOS,len(1)-start(1));
    for peptide=1:1:Row
        SelSeq{peptide}=SeqList{peptide}(start(peptide):len(peptide));
        indx_M=strfind(SelSeq{peptide},'M');AAs(peptide,1,indx_M)=1;
        indx_A=strfind(SelSeq{peptide},'A');AAs(peptide,2,indx_A)=1;
        indx_V=strfind(SelSeq{peptide},'V');AAs(peptide,3,indx_V)=1;
        indx_L=strfind(SelSeq{peptide},'L');AAs(peptide,4,indx_L)=1;
        indx_I=strfind(SelSeq{peptide},'I');AAs(peptide,5,indx_I)=1;
        indx_P=strfind(SelSeq{peptide},'P');AAs(peptide,6,indx_P)=1;
        indx_F=strfind(SelSeq{peptide},'F');AAs(peptide,7,indx_F)=1;
        indx_W=strfind(SelSeq{peptide},'W');AAs(peptide,8,indx_W)=1;
        indx_G=strfind(SelSeq{peptide},'G');AAs(peptide,9,indx_G)=1;
        indx_S=strfind(SelSeq{peptide},'S');AAs(peptide,10,indx_S)=1;
        indx_C=strfind(SelSeq{peptide},'C');AAs(peptide,11,indx_C)=1;
        indx_N=strfind(SelSeq{peptide},'N');AAs(peptide,12,indx_N)=1;
        indx_Q=strfind(SelSeq{peptide},'Q');AAs(peptide,13,indx_Q)=1;
        indx_Y=strfind(SelSeq{peptide},'Y');AAs(peptide,14,indx_Y)=1;
        indx_T=strfind(SelSeq{peptide},'T');AAs(peptide,15,indx_T)=1;
        indx_K=strfind(SelSeq{peptide},'K');AAs(peptide,16,indx_K)=1;
        indx_R=strfind(SelSeq{peptide},'R');AAs(peptide,17,indx_R)=1;
        indx_H=strfind(SelSeq{peptide},'H');AAs(peptide,18,indx_H)=1;
        indx_D=strfind(SelSeq{peptide},'D');AAs(peptide,19,indx_D)=1;
        indx_E=strfind(SelSeq{peptide},'E');AAs(peptide,20,indx_E)=1;
    end
elseif(strcmp(type,'int'))
    AAs=zeros(Row,len(1)-start(1));
    for peptide=1:1:Row
        SelSeq{peptide}=SeqList{peptide}(start(peptide):len(peptide));
        indx_M=strfind(SelSeq{peptide},'M');AAs(peptide,indx_M)=1;
        indx_A=strfind(SelSeq{peptide},'A');AAs(peptide,indx_A)=2;
        indx_V=strfind(SelSeq{peptide},'V');AAs(peptide,indx_V)=3;
        indx_L=strfind(SelSeq{peptide},'L');AAs(peptide,indx_L)=4;
        indx_I=strfind(SelSeq{peptide},'I');AAs(peptide,indx_I)=5;
        indx_P=strfind(SelSeq{peptide},'P');AAs(peptide,indx_P)=6;
        indx_F=strfind(SelSeq{peptide},'F');AAs(peptide,indx_F)=7;
        indx_W=strfind(SelSeq{peptide},'W');AAs(peptide,indx_W)=8;
        indx_G=strfind(SelSeq{peptide},'G');AAs(peptide,indx_G)=9;
        indx_S=strfind(SelSeq{peptide},'S');AAs(peptide,indx_S)=10;
        indx_C=strfind(SelSeq{peptide},'C');AAs(peptide,indx_C)=11;
        indx_N=strfind(SelSeq{peptide},'N');AAs(peptide,indx_N)=12;
        indx_Q=strfind(SelSeq{peptide},'Q');AAs(peptide,indx_Q)=13;
        indx_Y=strfind(SelSeq{peptide},'Y');AAs(peptide,indx_Y)=14;
        indx_T=strfind(SelSeq{peptide},'T');AAs(peptide,indx_T)=15;
        indx_K=strfind(SelSeq{peptide},'K');AAs(peptide,indx_K)=16;
        indx_R=strfind(SelSeq{peptide},'R');AAs(peptide,indx_R)=17;
        indx_H=strfind(SelSeq{peptide},'H');AAs(peptide,indx_H)=18;
        indx_D=strfind(SelSeq{peptide},'D');AAs(peptide,indx_D)=19;
        indx_E=strfind(SelSeq{peptide},'E');AAs(peptide,indx_E)=20;
    end
elseif(strcmp(type,'idx'))
    AAs=zeros(Row,len(1)-start(1));
    for peptide=1:1:Row
        SelSeq{peptide}=SeqList{peptide}(start(peptide):len(peptide));
        indx_M=strfind(SelSeq{peptide},'M');indx_M=(indx_M-1)*20+1;
        indx_A=strfind(SelSeq{peptide},'A');indx_A=(indx_A-1)*20+2;
        indx_V=strfind(SelSeq{peptide},'V');indx_V=(indx_V-1)*20+3;
        indx_L=strfind(SelSeq{peptide},'L');indx_L=(indx_L-1)*20+4;
        indx_I=strfind(SelSeq{peptide},'I');indx_I=(indx_I-1)*20+5;
        indx_P=strfind(SelSeq{peptide},'P');indx_P=(indx_P-1)*20+6;
        indx_F=strfind(SelSeq{peptide},'F');indx_F=(indx_F-1)*20+7;
        indx_W=strfind(SelSeq{peptide},'W');indx_W=(indx_W-1)*20+8;
        indx_G=strfind(SelSeq{peptide},'G');indx_G=(indx_G-1)*20+9;
        indx_S=strfind(SelSeq{peptide},'S');indx_S=(indx_S-1)*20+10;
        indx_C=strfind(SelSeq{peptide},'C');indx_C=(indx_C-1)*20+11;
        indx_N=strfind(SelSeq{peptide},'N');indx_N=(indx_N-1)*20+12;
        indx_Q=strfind(SelSeq{peptide},'Q');indx_Q=(indx_Q-1)*20+13;
        indx_Y=strfind(SelSeq{peptide},'Y');indx_Y=(indx_Y-1)*20+14;
        indx_T=strfind(SelSeq{peptide},'T');indx_T=(indx_T-1)*20+15;
        indx_K=strfind(SelSeq{peptide},'K');indx_K=(indx_K-1)*20+16;
        indx_R=strfind(SelSeq{peptide},'R');indx_R=(indx_R-1)*20+17;
        indx_H=strfind(SelSeq{peptide},'H');indx_H=(indx_H-1)*20+18;
        indx_D=strfind(SelSeq{peptide},'D');indx_D=(indx_D-1)*20+19;
        indx_E=strfind(SelSeq{peptide},'E');indx_E=(indx_E-1)*20+20;
        AAs=[indx_M indx_A indx_V indx_L indx_I indx_P indx_F indx_W indx_G indx_S indx_C indx_N indx_Q indx_Y indx_T indx_K indx_R indx_H indx_D indx_E];
    end

end
end





