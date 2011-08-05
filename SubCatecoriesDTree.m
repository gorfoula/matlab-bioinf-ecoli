function []=SubCatecoriesDTree(DataINfile,testFile)

fid=fopen(DataINfile,'r');
[PARAM] = textscan(fid,'%d',3,'delimiter','\t');
PEPTIDES_NUM_train=PARAM{1}(1);
FEATURES_NUM_train=PARAM{1}(2);

data=zeros(PEPTIDES_NUM_train,FEATURES_NUM_train);
read=cell(PEPTIDES_NUM_train,1);
for i=1:1:PEPTIDES_NUM_train
[read(i)] = textscan(fid,'%d',FEATURES_NUM_train,'delimiter','\t');
data(i,:)=read{i}(1:end);
end


fid_test=fopen(testFile,'r');
[PARAM] = textscan(fid_test,'%d',3,'delimiter','\t');
PEPTIDES_NUM=PARAM{1}(1);
FEATURES_NUM=PARAM{1}(2);

testdata=zeros(PEPTIDES_NUM,FEATURES_NUM);
testread=cell(PEPTIDES_NUM,1);
for i=1:1:PEPTIDES_NUM
[testread(i)] = textscan(fid_test,'%d',FEATURES_NUM,'delimiter','\t');
testdata(i,:)=testread{i}(1:end);
end
VALUES=20;
FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};
FeatureID=1:1:2000;
NumericCode=mod(FeatureID,VALUES);
PosOnPepseq=floor(FeatureID/VALUES)+(NumericCode>0);
NumericCode(NumericCode==0)=VALUES;
OneLetterCode=FEATURES(NumericCode);
NUMofFEATURES=FEATURES_NUM;
CYTO=double(sum(data(:,1)==0))./double(PEPTIDES_NUM_train);
PERI=double(sum(data(:,1)==1))./double(PEPTIDES_NUM_train);
OMBS=double(sum(data(:,1)==2))./double(PEPTIDES_NUM_train);
LIPO=double(sum(data(:,1)==3))./double(PEPTIDES_NUM_train);
%%%%%%%%%%%%%%%  Decision  Tree   %%%%%%%%%%%%%%%%%%%

Predictor_var=cell(NUMofFEATURES-1,1);
Canonical=FeatureID;
for i=1:1:NUMofFEATURES-1
    Predictor_var{i}=[OneLetterCode{i}(1:end),int2str(PosOnPepseq(i))];
end


Category_Names={'Cyto  ','Peri  ','OMBs  ','Lipo  '};
catg=Category_Names(data(:,1)+1)';
testcatg=Category_Names(testdata(:,1)+1)';
classprobs=[CYTO PERI OMBS LIPO];
% t = classregtree(data(:,2:end),catg,'names',Predictor_var,'prune','on','categorical',Canonical);

t = treefit(data(:,2:end),catg,'names',Predictor_var,'prune','on','catidx',Canonical,'method','classification','splitcriterion','gdi','priorprob',classprobs);

% [c,s,n,bestt] = treetest(t,'resubstitution');
[ct,st,nt,bestt] = treetest(t,'cross',testdata(:,2:end),testcatg,'nsamples',20);
if(bestt>0)
    tmin = prune(t,'level',bestt-1);
    treedisp(tmin,'names',Predictor_var);
    1-ct(bestt)
else
    tmin=t;
end
treedisp(t,'names',Predictor_var);



sfit = treeval(t,testdata(:,2:end));     % Find assigned class numbers
t_test=1-mean(abs(sfit-(testdata(:,1)+1)))  % Proportion in correct class

sfit_tmin = treeval(tmin,testdata(:,2:end));     % Find assigned class numbers
tmin_test=1-mean(abs(sfit_tmin-(testdata(:,1)+1)))  % Proportion in correct class

sfit = treeval(t,data(:,2:end));     % Find assigned class numbers
t_train=1-mean(abs(sfit-(data(:,1)+1)))  % Proportion in correct class

sfit_tmin = treeval(tmin,data(:,2:end));     % Find assigned class numbers
tmin_train=1-mean(abs(sfit_tmin-(data(:,1)+1)))  % Proportion in correct

end