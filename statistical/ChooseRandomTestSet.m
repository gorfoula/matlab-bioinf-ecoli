function []=ChooseRandomTestSet(IFile,perc)

found_back=regexp(IFile,'[\/]');
found_dot=regexp(IFile,'[.]');

dir=IFile(1:found_back(end));
filename=IFile(found_back(end)+1:found_dot(end)-1);

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file

proteins=length(GnNames);
testSetLen=ceil(proteins*perc);
rand_index=randperm(proteins);
test_index=rand_index(1:testSetLen);
train_index=rand_index(testSetLen+1:end);

WriteSelection([dir,'TEST/',filename,'_test.txt'],GnNames(test_index), Descr(test_index), TotLen(test_index), MW(test_index), SPLen(test_index), SEQ(test_index))
WriteSelection([dir,'TRAIN/',filename,'.txt'],GnNames(train_index), Descr(train_index), TotLen(train_index), MW(train_index), SPLen(train_index), SEQ(train_index))

end


function []=WriteSelection(File,GnNames, Descr, TotLen, MW, SPLen, SEQ)
s=fopen(File,'w');
for i=1:1:length(GnNames)   %for all best signal peptides
%     GnNames{i}
    invalid=regexp(SEQ{i},'[^MAVLIPFWGSCNQYTKRHDE]');
    if(isempty(invalid)==0)
        SEQ{i}(invalid)='A';
    end
    fprintf(s,'%s \t %s \t %d \t %d \t %d \t %s \n',GnNames{i}, Descr{i}, TotLen(i), ceil(molweight(SEQ{i})), SPLen(i), SEQ{i});  %% file save best SP in an index       
end
fclose all;

end
% 
% perc=0.3;
% ChooseRandomTestSet('Dataset/AllCytoplasmic_mn.txt',perc);
% ChooseRandomTestSet('Dataset/AllSecreted_mn.txt',perc);
% ChooseRandomTestSet('Dataset/Allb_barrel_mn.txt',perc);
% ChooseRandomTestSet('Dataset/AllPeriplasmic_mn.txt',perc);
% ChooseRandomTestSet('Dataset/AllLipoproteins_mn.txt',perc);