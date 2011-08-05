function [PhoMatrix Pho M]=Hydrophobicity(PeptideSeq,scale)
%% M: max value of the scale
%% Load Scale table
 [aa name Kyte Hopp Cornette Eisenberg Rose Janin Engelman Levitt WW1 WW2 WWd Phobic Philic] = textread('Features/PhoScale.txt','%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',-1,'delimiter','\t');
%% Assign Scale values acoording to scale choice
switch (scale)
    case 'K'
        Kyte=CellTable2Double(Kyte);
        PhoMatrix=(Kyte(PeptideSeq+1));
        M=max(abs(Kyte));
    case 'H'
        Hopp=CellTable2Double(Hopp);
        PhoMatrix=(Hopp(PeptideSeq+1));
        M=max(abs(Hopp));
    case 'C'
        Cornette=CellTable2Double(Cornette);
        PhoMatrix=(Cornette(PeptideSeq+1));
        M=max(abs(Cornette));
    case 'E'
        Eisenberg=CellTable2Double(Eisenberg);
        PhoMatrix=(Eisenberg(PeptideSeq+1));
        M=max(abs(Eisenberg));
    case 'R'
        Rose=CellTable2Double(Rose);
        PhoMatrix=(Rose(PeptideSeq+1));
        M=max(abs(Rose));
    case 'J'
        Janin=CellTable2Double(Janin);
        PhoMatrix=(Janin(PeptideSeq+1));
        M=max(abs(Janin));
    case 'En'
        Engelman=CellTable2Double(Engelman);
        PhoMatrix=(Engelman(PeptideSeq+1));
        M=max(abs(Engelman));
    case 'L'
        Levitt=CellTable2Double(Levitt);
        PhoMatrix=(Levitt(PeptideSeq+1));
        M=max(abs(Levitt));
    case 'WW1'
        WW1=CellTable2Double(WW1);
        PhoMatrix=(WW1(PeptideSeq+1));
        M=max(abs(WW1));
    case 'WW2'
        WW2=CellTable2Double(WW2);
        PhoMatrix=(WW2(PeptideSeq+1));
        M=max(abs(WW2));
    case 'WWd'
        WWd=CellTable2Double(WWd);
        PhoMatrix=(WWd(PeptideSeq+1));
        M=max(abs(WWd));
    case 'Phobic'
        Phobic=CellTable2Double(Phobic);
        PhoMatrix=(Phobic(PeptideSeq+1));
        M=max(abs(Phobic));
    case 'Philic'
        Philic=CellTable2Double(Philic);
        PhoMatrix=(Philic(PeptideSeq+1));    
        M=max(abs(Philic));
    otherwise
        display('!! Wrong scale name');
        PhoMatrix=PeptideSeq;
        return;
end

Pho=sum(PhoMatrix);
% PhoMatrix=PhoMatrix./sqrt((PhoMatrix.^2)+1);
PhoMatrix=PhoMatrix';
end