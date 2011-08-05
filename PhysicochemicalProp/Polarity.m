function [PolMatrix Pol M]=Polarity(PeptideSeq,scale)
%% M: max value of the scale
%AMINOS=20;
[aa name Zimmerman Grantham Levitt NetCharge NegCharge PosCharge Phobic] = textread('Features/PolScale.txt','%s %s %s %s %s %s %s %s %s',-1,'delimiter','\t');

switch (scale)
    case 'Z'
        Zimmerman=CellTable2Double(Zimmerman);
        PolMatrix=(Zimmerman(PeptideSeq+1));
        M=max(abs(Zimmerman));
    case 'G'
        Grantham=CellTable2Double(Grantham);
        PolMatrix=(Grantham(PeptideSeq+1));
        M=max(abs(Grantham));
    case 'L' %% alternative of Zimmer scale where less than 3 is negative
        Levitt=CellTable2Double(Levitt);
        PolMatrix=(Levitt(PeptideSeq+1));
        M=max(abs(Levitt));        
    case 'NC'
        NetCharge=CellTable2Double(NetCharge);
        PolMatrix=(NetCharge(PeptideSeq+1));
        M=max(abs(NetCharge));
    case 'Neg'
        NegCharge=CellTable2Double(NegCharge);
        PolMatrix=(NegCharge(PeptideSeq+1));
        M=max(abs(NegCharge));
    case 'Pos'
        PosCharge=CellTable2Double(PosCharge);
        PolMatrix=(PosCharge(PeptideSeq+1));
        M=max(abs(PosCharge));
    case 'Phobic'
        Phobic=CellTable2Double(Phobic);
        PolMatrix=(Phobic(PeptideSeq+1));
        M=max(abs(Phobic));
    otherwise
        display('!! Wrong scale name');
        PolMatrix=PeptideSeq;
        return;
end

Pol=sum(PolMatrix);
% PolMatrix=PolMatrix./sqrt((PolMatrix.^2)+1);
PolMatrix=PolMatrix';

end