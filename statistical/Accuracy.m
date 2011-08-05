function [ACC]=Accuracy(groups,allscores)

cp = classperf(double(groups));
Peptides=length(groups);
classes=ones(Peptides,1);
classes((allscores)<0)=0;
classperf(cp,classes,[1:1:Peptides]');
ACC=cp.CorrectRate;
end