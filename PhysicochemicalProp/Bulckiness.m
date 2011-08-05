function [BulkMatrix BulkSum]=Bulckiness(PeptideSeq)

%AMINOS=20;
[aa name Zimmerman] = textread('Features/Bulkiness.txt','%s %s %s',-1,'delimiter','\t');

BulkMatrix=str2double(Zimmerman(PeptideSeq+1))./max(str2double(Zimmerman));

BulkSum=sum(BulkMatrix);
BulkMatrix=BulkMatrix';

end