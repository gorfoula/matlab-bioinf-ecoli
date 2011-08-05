function []=SelectedFeatures(inputFILE,ftsel,DataINfile,testFile,train,polyorder)

%%%%%%%%%%%%%%%%%%  READ_FILES   %%%%%%%%%%%%%%%%%%%%%
[aa, FeatureID] = textread(inputFILE,'%d %d',-1,'delimiter','\t');
% aa=1:2000;
% FeatureID=1:2001;

NUMofFEATURES=aa(end);display(NUMofFEATURES);
FeatureID=[1;FeatureID];

% [data PEPTIDES_NUM_train]=ReadSelectedFeatures(DataINfile,aa,FeatureID);

fid=fopen(DataINfile,'r');
[PARAM] = textscan(fid,'%d',3,'delimiter','\t');
PEPTIDES_NUM_train=PARAM{1}(1);
FEATURES_NUM_train=PARAM{1}(2);
% FeatureID=[1;FeatureID];



data=zeros(PEPTIDES_NUM_train,aa(end)+1);
read=cell(PEPTIDES_NUM_train,1);
for i=1:1:PEPTIDES_NUM_train
[read(i)] = textscan(fid,'%d',FEATURES_NUM_train,'delimiter','\t');
data(i,:)=read{i}(FeatureID);
end

fid_test=fopen(testFile,'r');
[PARAM] = textscan(fid_test,'%d',3,'delimiter','\t');
PEPTIDES_NUM=PARAM{1}(1);
FEATURES_NUM=PARAM{1}(2);

testdata=zeros(PEPTIDES_NUM,aa(end)+1);
testread=cell(PEPTIDES_NUM,1);
for i=1:1:PEPTIDES_NUM
[testread(i)] = textscan(fid_test,'%d',FEATURES_NUM,'delimiter','\t');
testdata(i,:)=testread{i}(FeatureID);
end

fclose all;

%%%%%%%%%%%%%  TYPE OF REPRESSENTATION  %%%%%%%%%%%%%%
FeatureID=FeatureID(2:end);
FeatureID=FeatureID-1;
NumericCode=ones(NUMofFEATURES,1);
OneLetterCode=zeros(NUMofFEATURES,1);
PosOnPepseq=zeros(NUMofFEATURES,1);
OneLetterCode=cell(NUMofFEATURES,1);

switch (ftsel)
    case 0
        VALUES=20;
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};
    case 3
        VALUES=20;
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'Bulk' 'Pol' 'Pho'};
        index=(FeatureID<=3);
        index2=(FeatureID>3);
        FeatureID(index2)=FeatureID(index2)-3;
        NumericCode=mod(FeatureID(index2),VALUES);
        PosOnPepseq=floor(FeatureID(index2)/VALUES)+(NumericCode>0);
        NumericCode(NumericCode==0)=VALUES;
        
        FeatureID(index2)=FeatureID(index2)+(PosOnPepseq-1)*3;
        FeatureID(index)=FeatureID(index)+20;
        VALUES=23;
    case 4
        VALUES=20;
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'Bulk' 'PosCh' 'NegCh' 'Pho'};
        index=(FeatureID<=4);
        index2=(FeatureID>4);
        FeatureID(index2)=FeatureID(index2)-4;
        NumericCode=mod(FeatureID(index2),VALUES);
        PosOnPepseq=floor(FeatureID(index2)/VALUES)+(NumericCode>0);
        NumericCode(NumericCode==0)=VALUES;
        
        FeatureID(index2)=FeatureID(index2)+(PosOnPepseq-1)*4;
        FeatureID(index)=FeatureID(index)+20;
        VALUES=24;
    case 5
        VALUES=20;
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'SPLen' 'MATLen' 'Bulk' 'Pol' 'Pho'};
        index=(FeatureID<=5);
        index2=(FeatureID>5);
        FeatureID(index2)=FeatureID(index2)-5;
        NumericCode=mod(FeatureID(index2),VALUES);
        PosOnPepseq=floor(FeatureID(index2)/VALUES)+(NumericCode>0);
        NumericCode(NumericCode==0)=VALUES;
        
        FeatureID(index2)=FeatureID(index2)+(PosOnPepseq-1)*5;
        FeatureID(index)=FeatureID(index)+20;
        VALUES=25;
    case 6
        VALUES=20;
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'SPLen' 'MATLen' 'Bulk' 'PosCh' 'NegCh' 'Pho'};
        index=(FeatureID<=6);
        index2=(FeatureID>6);
        FeatureID(index2)=FeatureID(index2)-6;
        NumericCode=mod(FeatureID(index2),VALUES);
        PosOnPepseq=floor(FeatureID(index2)/VALUES)+(NumericCode>0);
        NumericCode(NumericCode==0)=VALUES;
        
        FeatureID(index2)=FeatureID(index2)+(PosOnPepseq-1)*6;
        FeatureID(index)=FeatureID(index)+20;
        VALUES=26;
    case 20119
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'ph' 'sml' 'R' 'K' 'D' 'E' 'x' 'pol' 'H' 'M' 'W' '@' '+' 'sm' 'o' 'h' 'b' 'P' 'q' 'C' 'Bulk' 'PosCh' 'NegCh' 'Pho'};
                
        index_prop=(FeatureID<=4);        
        index_20=and(FeatureID>4,FeatureID<=2004);
        index_11=and(FeatureID>2004,FeatureID<=3104);
        index_9=(FeatureID>3104);
        
        FeatureID(index_20)=FeatureID(index_20)-4;
        FeatureID(index_11)=FeatureID(index_11)-2004;
        FeatureID(index_9)=FeatureID(index_9)-3104;
               
        NumericCode(index_20)=mod(FeatureID(index_20),20);
        PosOnPepseq(index_20)=floor(FeatureID(index_20)/20)+(NumericCode(index_20)>0);
        NumericCode(and(NumericCode==0,index_20))=20;
        OneLetterCode(index_20)=FEATURES(NumericCode(index_20));
                      
        NumericCode(index_11)=mod(FeatureID(index_11),11);
        PosOnPepseq(index_11)=floor(FeatureID(index_11)/11)+(NumericCode(index_11)>0);
        NumericCode(and(NumericCode==0,index_11))=11;
        NumericCode(index_11)=NumericCode(index_11)+20;
        OneLetterCode(index_11)=FEATURES(NumericCode(index_11));
        
        NumericCode(index_9)=mod(FeatureID(index_9),9);
        PosOnPepseq(index_9)=floor(FeatureID(index_9)/9)+(NumericCode(index_9)>0);
        NumericCode(and(NumericCode==0,index_9))=9;
        NumericCode(index_9)=NumericCode(index_9)+31;
        OneLetterCode(index_9)=FEATURES(NumericCode(index_9));
        
        NumericCode(index_prop)=FeatureID(index_prop)+40;
        PosOnPepseq(index_prop)=-1;
        OneLetterCode(index_prop)=FEATURES(NumericCode(index_prop));
    case 201192
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'ph' 'sml' 'R' 'K' 'D' 'E' 'x' 'pol' 'H' 'M' 'W' '@' '+' 'sm' 'o' 'h' 'b' 'P' 'q' 'C'};
        
        index_20=(FeatureID<=2000);
        index_11=and(FeatureID>2000,FeatureID<=3100);
        index_9=(FeatureID>3100);
        
        FeatureID(index_11)=FeatureID(index_11)-2000;
        FeatureID(index_9)=FeatureID(index_9)-3100;
               
        NumericCode(index_20)=mod(FeatureID(index_20),20);
        PosOnPepseq(index_20)=floor(FeatureID(index_20)/20)+(NumericCode(index_20)>0);
        NumericCode(and(NumericCode==0,index_20))=20;
        OneLetterCode(index_20)=FEATURES(NumericCode(index_20));
                      
        NumericCode(index_11)=mod(FeatureID(index_11),11);
        PosOnPepseq(index_11)=floor(FeatureID(index_11)/11)+(NumericCode(index_11)>0);
        NumericCode(and(NumericCode==0,index_11))=11;
        NumericCode(index_11)=NumericCode(index_11)+20;
        OneLetterCode(index_11)=FEATURES(NumericCode(index_11));
        
        NumericCode(index_9)=mod(FeatureID(index_9),9);
        PosOnPepseq(index_9)=floor(FeatureID(index_9)/9)+(NumericCode(index_9)>0);
        NumericCode(and(NumericCode==0,index_9))=9;
        NumericCode(index_9)=NumericCode(index_9)+31;
        OneLetterCode(index_9)=FEATURES(NumericCode(index_9));
    case 201193
        FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'ph' 'sml' 'R' 'K' 'D' 'E' 'x' 'pol' 'H' 'M' 'W' '@' '+' 'sm' 'o' 'h' 'b' 'P' 'q' 'C' 'Pos' 'Pho'};
        
        index_prop=(FeatureID<=2);
        index_20=and(FeatureID>2,FeatureID<=2002);
        index_11=and(FeatureID>2002,FeatureID<=3102);
        index_9=(FeatureID>3102);
        
        FeatureID(index_20)=FeatureID(index_20)-2;
        FeatureID(index_11)=FeatureID(index_11)-2002;
        FeatureID(index_9)=FeatureID(index_9)-3102;
               
        NumericCode(index_20)=mod(FeatureID(index_20),20);
        PosOnPepseq(index_20)=floor(FeatureID(index_20)/20)+(NumericCode(index_20)>0);
        NumericCode(and(NumericCode==0,index_20))=20;
        OneLetterCode(index_20)=FEATURES(NumericCode(index_20));
                      
        NumericCode(index_11)=mod(FeatureID(index_11),11);
        PosOnPepseq(index_11)=floor(FeatureID(index_11)/11)+(NumericCode(index_11)>0);
        NumericCode(and(NumericCode==0,index_11))=11;
        NumericCode(index_11)=NumericCode(index_11)+20;
        OneLetterCode(index_11)=FEATURES(NumericCode(index_11));
        
        NumericCode(index_9)=mod(FeatureID(index_9),9);
        PosOnPepseq(index_9)=floor(FeatureID(index_9)/9)+(NumericCode(index_9)>0);
        NumericCode(and(NumericCode==0,index_9))=9;
        NumericCode(index_9)=NumericCode(index_9)+31;
        OneLetterCode(index_9)=FEATURES(NumericCode(index_9));
        
        NumericCode(index_prop)=FeatureID(index_prop)+40;
        PosOnPepseq(index_prop)=-1;
        OneLetterCode(index_prop)=FEATURES(NumericCode(index_prop));                             
    case 9
        VALUES=9;
        FEATURES={'@' '+' 's' 'o' 'h' 'b' 'P' 'q' 'C'};
    case 10
        VALUES=10;
        FEATURES={'@' '+' 's' 'o' 'h' 'b' 'P' 'q' 'C' 'A'};
    case 11
        VALUES=11;
        FEATURES={'h' 's' 'R' 'K' 'D' 'E' 'x' 'q' 'H' 'M' 'W'}; 
%         FEATURES={'pho' 'small' 'R' 'K' 'D' 'E' 'Hyx' 'q' 'H' 'M' 'W'}; 
    otherwise
        display('=>Wrong type of characteristics at <SelectedFeatures.m>');
        return;
end
%%%%%%%%%%%%%%  CORRESPONDING FEATURES  %%%%%%%%%%%%%%%
if(ftsel<20)
NumericCode=mod(FeatureID,VALUES);
PosOnPepseq=floor(FeatureID/VALUES)+(NumericCode>0);
NumericCode(NumericCode==0)=VALUES;
OneLetterCode=FEATURES(NumericCode);
PosOnPepseq(NumericCode>20)=-PosOnPepseq(NumericCode>20);
end

% highef=data(data(:,1)==0,:);
% medef=data(data(:,1)==1,:);
% lowef=data(data(:,1)==2,:);
% 
% data=[lowef;highef];

%%%%%%%%%%%%%%%%  STATISTICAL FEATURES %%%%%%%%%%%%%
class1=(data(:,1)==max(data(:,1)));
class2=(data(:,1)==min(data(:,1)));
Periplasmic=data(class1,2:end);
Cytoplasmic=data(class2,2:end);

MeanAll=mean(data(:,2:end));

PeriplasmicFeatureMean=mean(Periplasmic);
CytoplasmicFeatureMean=mean(Cytoplasmic);

derak=zeros(PEPTIDES_NUM_train,NUMofFEATURES);
derak(:,1)=1;
all_max=conv2(derak,max(abs(data(:,2:end))));
all_max=all_max(:,1:NUMofFEATURES);

PeriplasmicFeatureMean_norm=mean(Periplasmic./all_max(class1,:));
CytoplasmicFeatureMean_norm=mean(Cytoplasmic./all_max(class2,:));

CorCoef=tril(corr(double(data(:,2:end))),-1);
[x,y]=find(CorCoef>=0.5);
coordinates=[x y CorCoef(CorCoef>=0.5) MeanAll(x)' MeanAll(y)' PosOnPepseq(x) PosOnPepseq(y)];
AminoCouple=[OneLetterCode(x) OneLetterCode(y)]  %%%

PeriplasmicCorCoef=tril(corr(double(Periplasmic)),-1);
CytoplasmicCorCoef=tril(corr(double(Cytoplasmic)),-1);

[x,y]=find(PeriplasmicCorCoef>=0.5);
coordinates=[x y PeriplasmicCorCoef(PeriplasmicCorCoef>=0.5) PeriplasmicFeatureMean(x)' PeriplasmicFeatureMean(y)' PosOnPepseq(x) PosOnPepseq(y)];
AminoCouple=[OneLetterCode(x)' OneLetterCode(y)']

%%%%%%%%%%%%%%  Decision  Tree   %%%%%%%%%%%%%%%%%%%
% POS=Double2CellTable(PosOnPepseq);
% [Features_Names]=MergeColumns([OneLetterCode' POS],'_');
% header=MergeColumns([{'Label'} Features_Names'],';');
% 
% data_cell=Double2CellTable(testdata(:,2:end));
% body=MergeColumns(data_cell,';');
% 
% Category_Names={'Cyto','Secr'};
% catg=Category_Names(testdata(:,1)+1)';
% body=MergeColumns([catg body],';');
% FileWriteTable('Gems/Curated/TEST/CytoVSsecr_mat_sel_test.csv',body,header{1,1},'w');
% 
% 
% Predictor_var=zeros(NUMofFEATURES,4);
% temp=int2str(PosOnPepseq);
% Canonical=aa(NumericCode<21);
% for i=1:1:NUMofFEATURES
%     if(NumericCode(i)<21)
%         Predictor_var{i}=[OneLetterCode{i}(1:end),temp(i,:)];
%     else
%         Predictor_var{i}=OneLetterCode{i}(1:end);
%     end
% end
% 
% Category_Names={'Cyto','Secr'};
% catg=Category_Names(testdata(:,1)+1)';
% 
% 
% testcatg=Category_Names(testdata(:,1)+1)';
% 
% t = classregtree(data(:,2:end),catg,'names',Predictor_var,'prune','on','categorical',Canonical);
% 
% t = treefit(data(:,2:end),catg,'names',Predictor_var,'prune','on','catidx',Canonical,'method','classification','splitcriterion','gdi','cost',[0 1;1 0]);
% 
% [c,s,n,best] = treetest(t,'resubstitution');
% [ct,st,nt,bestt] = treetest(t,'cross',testdata(:,2:end),testcatg,'nsamples',20);
% if(bestt>0)
%     tmin = prune(t,'level',bestt-1);
%     treedisp(tmin,'names',Predictor_var);
%     1-ct(bestt)
% else
%     tmin=t;
% end
% treedisp(t,'names',Predictor_var);
% 
% 
% 
% sfit = treeval(t,testdata(:,2:end));     % Find assigned class numbers
% t_test=1-mean(abs(sfit-(testdata(:,1)+1)))  % Proportion in correct class
% 
% sfit_tmin = treeval(tmin,testdata(:,2:end));     % Find assigned class numbers
% tmin_test=1-mean(abs(sfit_tmin-(testdata(:,1)+1)))  % Proportion in correct class
% 
% sfit = treeval(t,data(:,2:end));     % Find assigned class numbers
% t_train=1-mean(abs(sfit-(data(:,1)+1)))  % Proportion in correct class
% 
% sfit_tmin = treeval(tmin,data(:,2:end));     % Find assigned class numbers
% tmin_train=1-mean(abs(sfit_tmin-(data(:,1)+1)))  % Proportion in correct class

%%%%%%%%%%%%%%%  SVM TRAINING   %%%%%%%%%%%%%%%%%%%
groups=data(:,1);
if(train==0)
    [w] = textread([inputFILE,'w.txt'],'%f',-1,'delimiter','\t');
else
    [train, test] = crossvalind('holdOut',groups);
    cp = classperf(groups);

    options = optimset('Display','iter','MaxIter',10^8);
    svmStruct = svmtrain(double(data(train,2:end)),groups(train),'Kernel_Function','polynomial','POLYORDER',polyorder,'METHOD','LS','QUADPROG_OPTS',options,'AUTOSCALE', false);
%     svmStruct = svmtrain(double(data(train,2:end)),groups(train),'Kernel_Function','polynomial','POLYORDER',polyorder,'METHOD','QP','QUADPROG_OPTS',options,'AUTOSCALE', false);
%     svmStruct = svmtrain(double(data(:,2:end)),data(:,1),'Kernel_Function','polynomial','POLYORDER',polyorder,'METHOD','QP','QUADPROG_OPTS',options,'AUTOSCALE', false);

%     opts = svmsmoset('Display','final','MaxIter',10^6);
%     svmStruct = svmtrain(double(data(:,2:end)),data(:,1),'Kernel_Function','polynomial','POLYORDER',polyorder,'METHOD','SMO','AUTOSCALE', false,'SMO_Opts',opts);
    [w]=getWeightFactors(svmStruct);
    dlmwrite([inputFILE,'w.txt'], w,'delimiter','\t');  %  write weights
%     classes = svmclassify(svmStruct,data(test,2:end));
%     classperf(cp,classes,test);
%     cp.CorrectRate;
%     AUC(classes, data(test,1), 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


groups=testdata(:,1);
cp = classperf(groups);
Peptides=length(groups);
score=testdata(:,2:end)*w(1:end-1)+w(end);
% for i=0.01:.01:.73
%     w_=[w(1:end-1)-i;w(end)];
%     score=testdata(:,2:end)*w_(1:end-1)+w_(end);
%     display(num2str(sum(score(groups==0)<0)));
%     display(num2str(sum(score(groups==0)>0)));
%     display(num2str(sum(score(groups==1)<0)));
%     display(num2str(sum(score(groups==1)>0)));
% end
classes=ones(Peptides,1);
classes(score<0)=0;
classperf(cp,classes,[1:1:Peptides]');
display(['Accuracy TEST: ',num2str(cp.CorrectRate)]);
allscores_=score+max(abs(score));
allscores=allscores_./max(allscores_);
display(['AUC TEST: ',num2str(AUC(groups, allscores, 1))]);

% FigureLegends(unique(score),score,groups,[{''} {'Hydrophobic Island Lenght'} {'Percent over all islands'}],[{'Peri'} {'Cyto'}],'p',{'-','.';':','o'},0);

groups=data(:,1);
cp = classperf(groups);
Peptides=length(groups);
score=data(:,2:end)*w(1:end-1)+w(end);
classes=ones(Peptides,1);
classes(score<0)=0;
classperf(cp,classes,[1:1:Peptides]');
display(['Accuracy TRAIN: ',num2str(cp.CorrectRate)]);
allscores_=score+max(abs(score));
allscores=allscores_./max(allscores_);
display(['AUC TRAIN: ',num2str(AUC(groups, allscores, 1))]);


%%%%%%%%%%%%   WRITE STATISTIC FILE   %%%%%%%%%%%%%

SFGraph(inputFILE,[w(1:end-1) NumericCode PosOnPepseq PeriplasmicFeatureMean_norm' CytoplasmicFeatureMean_norm' MeanAll'],FEATURES);


%%%%%%%%%%   WRITE CORRESPONDANCE FILE   %%%%%%%%%%%
s=fopen([inputFILE,'.cor.txt'],'w');
fprintf(s,'W\tFeatureID\tNumericCode\tPosOnPepseq\tOneLetterCode\tMeanOvePeriplasmic\tMeanOverCytoplasmic\tMeanOverAll\n');
for i=1:1:NUMofFEATURES
    fprintf(s,'%f\t%d\t%d\t%d\t%s\t',w(i),FeatureID(i)+1,NumericCode(i),PosOnPepseq(i),OneLetterCode{i});
    fprintf(s,'%1.3f\t%1.3f\t%1.3f\n',PeriplasmicFeatureMean(i),CytoplasmicFeatureMean(i),MeanAll(i));
end
fclose(s);

s=fopen([inputFILE,'.stc.txt'],'w');
fprintf(s,'FeatureID1\tFeatureID2\tCorCoef\tMeanFeature1\tMeanFeature2\tSeqPos1\tSeqPos2\tOneLettCode1\tOneLettCode2\n');

for i=1:1:length(x)
    fprintf(s,'%d\t%d\t%1.3f\t%1.3f\t%1.3f\t',coordinates(i,1),coordinates(i,2),coordinates(i,3),coordinates(i,4),coordinates(i,5));
    fprintf(s,'%d\t%d\t%s\t%s\n',coordinates(i,6),coordinates(i,7),AminoCouple{i,1},AminoCouple{i,2});
end
fclose all;


fclose all;
end