function [SCORES groups]=Prediction(DBIFile,ModelNames,TestSet,sel,gems,wSelFiles,outDir)
%% INPUT:
%%  ModelNames: SP and Mature Model Name (in order to read files read files)
%%  TestSet:    Pathway of Secreted and Cytoplasmic TEST sets
%%  gems:       1   for writing gems files
%%              2   for cytoplasmic file to be the RANDOM one
%%  sel:        if features are grouped (20119) or not (0)
%%  mat:        1   in order to take only mature domain of cytoplasmic
%%  wSelFiles:  1   write Secreted and Non Secreted proteins selected from
%%                  DBIFile
%%              2   Write files and split randomly into TEST and TRAIN set
%%  outDir:     directory of models and of output files if <wSelFiles> is 1 
close all;
%% Read Database File
[text_] = textread(DBIFile,'%s',-1,'delimiter','\n'); % read all lines
header_long=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(2:end),'[\t]'); % split coloumns
Table=CellTable2StrTable(Table_); % convert Cell table to string table
%% Read Selected Features File
[aa, FeatureID_sp] = textread([outDir,ModelNames{1},'.txt'],'%d %d',-1,'delimiter','\t');
FeatureID_sp=FeatureID_sp-1;  %%%% we don't want first column which is for category
[aa, FeatureID_mat] = textread([outDir,ModelNames{2},'.txt'],'%d %d',-1,'delimiter','\t');
FeatureID_mat=FeatureID_mat-1;  %%%% we don't want first column which is for category
%% Read weight file
[w_sp] = textread([outDir,ModelNames{1},'.txtw.txt'],'%f',-1,'delimiter','\t');  %%% read weights
[w_mat] = textread([outDir,ModelNames{2},'.txtw.txt'],'%f',-1,'delimiter','\t');  %%% read weights
VALUES=20;  %%% division factor / number of attribute IDs
len=100;
%% INIT OUTPUT TABLES (-700 is a score not given as an output by separation equation)
SCORES=ones(length(Table),2)*(-700);
GROUPS=ones(length(Table),1)*(-700);
%% DEFINITION OF SOME INDEXES
index_hm=CellTable2Double(Table(:,12))>0.7;
index_tat=(strcmp(Table(:,9),'TAT')) | (strcmp(Table(:,9),'TAT / SEC'));
index_Sec=(strcmp(Table(:,9),'SEC')) | (strcmp(Table(:,9),'TAT / SEC'));  %% select SEC system secreted proteins
index_cyto=strcmp(Table(:,8),'A') | strcmp(Table(:,8),'F1') | strcmp(Table(:,8),'A_trl');
index_lipo=(strcmp(Table(:,8),'I') | strcmp(Table(:,8),'E')) & strcmp(Table(:,9),'SEC');
index_bb=strcmp(Table(:,8),'H') & strcmp(Table(:,9),'SEC');
index_peri=strcmp(Table(:,8),'G') & strcmp(Table(:,9),'SEC');
index_peripheral=strcmp(Table(:,8),'F1');
index_Sec=and(index_Sec,not(strcmp(Table(:,8),'B')));  % no Inner Membrane (B) special case of proteins

TAT=Table(index_tat,:);
SP=Table(index_tat,23);
SEQ=Table(index_tat,25);
twinMotix=zeros(1,length(SP));
SPseq=cell(1,length(SP));
for i=1:size(TAT,1)
   curSP= SEQ{i}(1:str2double(SP{i}));
   [found]=regexp(curSP,'[R]{2}');
   twinMotix(i)=length(found);
   display([curSP,' twin ',num2str(length(found))]);
   SPseq{i}=curSP;
end

FileWriteTable('TATTable.txt',[TAT(:,1) SPseq' Double2CellTable(twinMotix)'],[],'w');

%% DEFINITION OF INDEXES ACOORDING TO TEST SET
% [INDEX_Sec SPLEN_Sec]=MatchNames(Table(:,1),TestSet{1});
% [INDEX_cyto SPLEN_cyto]=MatchNames(Table(:,1),TestSet{2});
%% Define cleavage
[SPLenL]=ChooseCleavage(Table(:,8),Table(:,9),Table(:,11),Table(:,13),Table(:,15),Table(:,23),Table(:,21));
% SPLenL(INDEX_Sec)=SPLEN_Sec(INDEX_Sec);
% SPLenL(INDEX_cyto)=SPLEN_cyto(INDEX_cyto);
% SPLenL=SPLenL';
%% Select proteins that have enough ammino acids for further analysis
TotLen=CellTable2Double(Table(:,5));

index_Len=(TotLen>SPLenL+len+1) & (TotLen>len);

INDEX=(index_lipo | index_cyto | index_bb | index_peri) & index_Len;
% INDEX=and(INDEX_Sec,index_Len);
GROUPS(index_cyto & index_Len )=0;
GROUPS((index_lipo | index_bb | index_peri) & index_Len)=1;

%% Select proteins that have enough ammino acids for further analysis
SPLEN=Double2CellTable(SPLenL(INDEX));
Names=Table(INDEX,1);   %Gene Names
DESCR=Table(INDEX,4);   %Description
LOCI=Table(INDEX,8);
HMM=Table(INDEX,12);
NN=Table(INDEX,14);
LIPO=Table(INDEX,17);
SecrSys=Table(INDEX,9);
SELECTED=Table(INDEX,end);    %Select Petides with Seq long enough
TOTLEN=Table(INDEX,5);
%% Define Secreted and Cytoplasmic
groups=logical(GROUPS(GROUPS>-700));
groups_secr=index_Sec & index_Len;groups_secr=groups_secr(GROUPS>-700);
groups_cyto=index_cyto & index_Len;groups_cyto=groups_cyto(GROUPS>-700);
groups_peri=index_peri & index_Len;groups_peri=groups_peri(GROUPS>-700);
groups_lipo=index_lipo & index_Len;groups_lipo=groups_lipo(GROUPS>-700);
groups_bb=index_bb & index_Len;groups_bb=groups_bb(GROUPS>-700);

Table_=Table(INDEX,:);
[allscores_mat]=ScoreCalculation(SELECTED,SPLenL(INDEX),w_mat,FeatureID_mat,VALUES,len,sel);   %%%% Score calculation
[allscores_sp]=ScoreCalculation(SELECTED,[],w_sp,FeatureID_sp,VALUES,len,sel);   %%%% Score calculation
SCORES=[allscores_sp allscores_mat];
%% AUC of our Models
display(['<Accurancy> SP: ',num2str(Accuracy(groups,allscores_sp)),'|  MAT: ',num2str(Accuracy(groups,allscores_mat))]);
display(['<AUC> SP: ',num2str(AUC(groups, NormalizeScores(allscores_sp), 1)),'| MAT: ',num2str(AUC(groups, NormalizeScores(allscores_mat), 1))]);
%% OTHER TOOLS
allscores_=(CellTable2Double(LIPO))./max((CellTable2Double(LIPO)));
display(['AUC LIPO: ',num2str(AUC(groups, allscores_, 1))]);
allscores_=(CellTable2Double(HMM))./max(CellTable2Double(HMM));
display(['AUC HMM: ',num2str(AUC(groups, allscores_, 1))]);
allscores_=(CellTable2Double(NN))./max(CellTable2Double(NN));
display(['AUC NN: ',num2str(AUC(groups, allscores_, 1))]);

header='GnNames\tDescr\tTotLen\tMW\tSPLen\tSEQ';
filenames={'AllCyto','AllSecr','AllLipo','AllPeri','Allbb','RandomSeq'};
legendnames={'Cytoplasmic','Secreted','Periplasmic','Lipoproteins','outer mem b-barrel','Cytoplasmic'};
legendnames_={'Lipoproteins','Periplasmic','Cytoplasmic','outer mem b-barrel'};
%% SPLEN and mature Distributions

%% Figure titles and legends
    SPLEN_=CellTable2Double(SPLEN);
    splen=min(SPLEN_):4:max(SPLEN_);
    [h1]=FigureLegends(splen,Freqcalc(SPLEN_(groups_secr),splen)',1,'Length (aas)','Percent','Signal Peptide length',legendnames{2},'b');
    display(['Mean SP: ',num2str(mean(SPLEN_(groups_secr)))]);
    
    lipo_sum=sum(groups_lipo);
    peri_sum=sum(groups_peri);
    bb_sum=sum(groups_bb);
    cyto_sum=sum(groups_cyto);
    secr_sum=sum(groups_secr);
    
    splen=12:4:max(SPLEN_);
    freq_sp=[Freqcalc(SPLEN_(groups_lipo),splen)' Freqcalc(SPLEN_(groups_peri),splen)' Freqcalc(SPLEN_(groups_bb),splen)'];
    truey=[SPLEN_(groups_lipo);SPLEN_(groups_peri);SPLEN_(groups_bb)];
    groups=[ones(1,lipo_sum) ones(1,peri_sum)*2  ones(1,bb_sum)*3];
    [h2]=FigureLegends(splen,freq_sp,2,'Length (aas)','Percent','Signal Peptide length',{'Lipoproteins','Periplasmic','outer mem b-barrel'},'b',{'-' '.'},truey,groups);
    display(['Mean SP LIPO : ',num2str(mean(SPLEN_(groups_lipo)))]);
    display(['Mean SP PERI : ',num2str(mean(SPLEN_(groups_peri)))]);
    display(['Mean SP BB : ',num2str(mean(SPLEN_(groups_bb)))]);

    
    MATLEN_=CellTable2Double(TOTLEN)-SPLEN_;
    matlen=min(MATLEN_):100:1000;
    freq_mat=[Freqcalc(MATLEN_(groups_lipo),matlen)' Freqcalc(MATLEN_(groups_peri),matlen)' Freqcalc(CellTable2Double(TOTLEN(groups_cyto)),matlen)' Freqcalc(MATLEN_(groups_bb),matlen)'];
    truey=[MATLEN_(groups_lipo);MATLEN_(groups_peri);MATLEN_(groups_cyto);MATLEN_(groups_bb)];
    groups=[ones(1,lipo_sum) ones(1,peri_sum)*2  ones(1,cyto_sum)*3 ones(1,bb_sum)*4];    
    [h3]=FigureLegends(matlen,freq_mat,3,'Length (aas)','Percent','Mature Domain length',legendnames_,'b',{'-'},truey,groups);
    
    freq_mat=[Freqcalc(CellTable2Double(TOTLEN(groups_cyto)),matlen)' Freqcalc(MATLEN_(groups_secr),matlen)'];
    truey=[CellTable2Double(TOTLEN(groups_cyto));MATLEN_(groups_secr)];
    groups=[ones(1,cyto_sum) ones(1,secr_sum)*2];
    [h4]=FigureLegends(matlen,freq_mat,4,'Length (aas)','Percent','Signal Peptide length',legendnames(1:2),'b',{'-'},truey,groups);



%% SELECTED FROM DATABASE /  WRITE IN FILES
if ((wSelFiles==1) || (wSelFiles==2))
    display('Writing Files of selected subcategories.');
%         FileWriteTable([outDir,filenames{1},'_alltools.txt'],Table_(groups,:),header_long,'w');
%         FileWriteTable([outDir,filenames{2},'_alltools.txt'],Table_(not(groups),:),header_long,'w');
    FileWriteTable([outDir,filenames{1},'.txt'],[Names(groups_cyto) DESCR(groups_cyto) TOTLEN(groups_cyto) LOCI(groups_cyto) SPLEN(groups_cyto) SELECTED(groups_cyto)],[],'w');
    FileWriteTable([outDir,filenames{2},'.txt'],[Names(groups_secr) DESCR(groups_secr) TOTLEN(groups_secr) LOCI(groups_secr) SPLEN(groups_secr) SELECTED(groups_secr)],[],'w');
    FileWriteTable([outDir,filenames{3},'.txt'],[Names(groups_lipo) DESCR(groups_lipo) TOTLEN(groups_lipo) LOCI(groups_lipo) SPLEN(groups_lipo) SELECTED(groups_lipo)],[],'w');
    FileWriteTable([outDir,filenames{4},'.txt'],[Names(groups_peri) DESCR(groups_peri) TOTLEN(groups_peri) LOCI(groups_peri) SPLEN(groups_peri) SELECTED(groups_peri)],[],'w');
    FileWriteTable([outDir,filenames{5},'.txt'],[Names(groups_bb) DESCR(groups_bb) TOTLEN(groups_bb) LOCI(groups_bb) SPLEN(groups_bb) SELECTED(groups_bb)],[],'w');
end
if(wSelFiles==2)
        ChooseRandomTestSet([outDir,filenames{1},'.txt'],.2);
        ChooseRandomTestSet([outDir,filenames{2},'.txt'],.2);
        ChooseRandomTestSet([outDir,filenames{3},'.txt'],.2);
        ChooseRandomTestSet([outDir,filenames{4},'.txt'],.2);
        ChooseRandomTestSet([outDir,filenames{5},'.txt'],.2);
end
%% GEMS OUTPUT %%
switch (gems)
    case 2
        FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TRAIN/',filenames{3},'.txt']} ; {[outDir,'TRAIN/',filenames{3},'.txt']}],[],'w');
        AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TRAIN/RandVSsecr_mat.gemsin'],'bin',100,1,0,0,1,[outDir,'FASTA/RandVSsecr_mat.fasta'],0);
        FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TEST/',filenames{3},'_test.txt']} ; {[outDir,'TEST/',filenames{3},'_test.txt']}],[],'w');
        AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TEST/RandVSsecr_mat.test'],'bin',100,1,0,0,1,[outDir,'FASTA/RandVSsecr_test_mat.fasta'],0);
    case 1

            FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TRAIN/',filenames{1},'.txt']} ; {[outDir,'TRAIN/',filenames{2},'.txt']}],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TRAIN/CytoVSsecr_S.gemsin'],'bin',52,-1,0,0,1,[outDir,'FASTA/CytoVSsecr_S.fasta'],0);
            FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TEST/',filenames{1},'_test.txt']} ; {[outDir,'TEST/',filenames{2},'_test.txt']}],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TEST/CytoVSsecr_S.test'],'bin',52,-1,0,0,1,[outDir,'FASTA/CytoVSsecr_S_test.fasta'],0);

            FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TRAIN/',filenames{1},'.txt']} ; {[outDir,'TRAIN/',filenames{2},'.txt']}],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TRAIN/CytoVSsecr_S_mat.gemsin'],'bin',100,1,0,0,1,[outDir,'FASTA/CytoVSsecr_S_mat.fasta'],0);
            FileWriteTable([outDir,'Input_files.txt'],[{[outDir,'TEST/',filenames{1},'_test.txt']} ; {[outDir,'TEST/',filenames{2},'_test.txt']}],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TEST/CytoVSsecr_S_mat.test'],'bin',100,1,0,0,1,[outDir,'FASTA/CytoVSsecr_S_mat_test.fasta'],0);

            FileWriteTable([outDir,'Input_files.txt'],[ {[outDir,'TEST/',filenames{1},'_test.txt']} ; {[outDir,'TEST/',filenames{2},'_test.txt']} ],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TEST/CytoVSsecr_S_20119.test'],'bin',100,-1,0,0,1,[outDir,'FASTA/CytoVSsecr_S_20119_test.fasta'],20119);
            FileWriteTable([outDir,'Input_files.txt'],[ {[outDir,'TRAIN/',filenames{1},'.txt']} ; {[outDir,'TRAIN/',filenames{2},'.txt']} ],[],'w');
            AllCategoriesCoding([outDir,'Input_files.txt'],[outDir,'TRAIN/CytoVSsecr_S_20119.gemsin'],'bin',100,-1,0,0,1,[outDir,'FASTA/CytoVSsecr_S_20119.fasta'],20119);            
    otherwise
        display('<gems> parameter 1 for input files 2 for cyto to be random');
end

% [Thres_hmm]=DefineScoreThres(CellTable2Double(Table(index_Len,10)),groups,groups_cyto,groups_peripheral,'SignalP HMM Smax Distribution',1);
% [Thres_nn]=DefineScoreThres(CellTable2Double(Table(index_Len,12)),groups,groups_cyto,groups_peripheral,'SignalP NN Dscore Distribution',2);
% [Thres_lipo]=DefineScoreThres(CellTable2Double(Table(index_Len,15)),groups,groups_cyto,groups_peripheral,'LipoP Score Distribution',3);
% [Thres_tmhhm]=DefineScoreThres(CellTable2Double(Table(index_Len,17)),groups,groups_cyto,groups_peripheral,'TMHMM aas in helix Distribution',4);
% [Thres_phob]=DefineScoreThres(CellTable2Double(Table(index_Len,20)),groups,groups_cyto,groups_peripheral,'Phobious No of helix Distribution',5);
% [Thres_tat]=DefineScoreThres(CellTable2Double(Table(index_Len,22)),groups,groups_cyto,groups_peripheral,'TatP Score Distribution',8);

% High_nn=and(or(CellTable2Double(Table(index_Len,12))>Thres_nn(2),CellTable2Double(Table(index_Len,17))>18),or(groups_peripheral,groups_cyto));
% 
% FileWriteTable([outDir,'Contradictional_alltools.txt'],Table_(High_nn,:),header_long,'w');
% FileWriteTable([dirData,'_OUT_IN_prediction.txt'],[Names(High_nn) Double2CellTable(allscores(High_nn))],'Name\tScore','w');
% [dirData,'_OUT_IN_prediction.txt']

end

function [SPLenL]=ChooseCleavage(LOCI,SecrSys,SignalP,SignalPNN,LipoP,TargetP,Phobious)

tokens=length(LOCI);
SPLenL=zeros(tokens,1);
%% Tat Cleavage
Tat=strcmp(SecrSys,'TAT / SEC');
SPLenL(Tat)=CellTable2Double(TargetP(Tat));
%% Sec Cleavage
% Sec=or(or(or(or(or(or(strcmp(LOCI,'E'),strcmp(LOCI,'I')),strcmp(LOCI,'G')),strcmp(LOCI,'H')),strcmp(LOCI,'X')),strcmp(LOCI,'F2')),strcmp(LOCI,'F3'));\
Sec=strcmp(SecrSys,'SEC');
sig=CellTable2Double(SignalP(Sec))-1;
signn=CellTable2Double(SignalPNN(Sec))-1;
phob=CellTable2Double(Phobious(Sec));
lip=CellTable2Double(LipoP(Sec));
sig(sig==0)=signn(sig==0);
sig(sig==0)=phob(sig==0);
sig(sig<=10)=lip(sig<=10);
SPLenL(Sec)=sig;
%% Lipo Cleavage
Lip=or(strcmp(LOCI,'I'),strcmp(LOCI,'E'));
SPLenL(Lip)=CellTable2Double(LipoP(Lip));

end

function [INDEX SPLEN]=MatchNames(DatabaseNames,IFile)

[text_] = textread(IFile,'%s',-1,'delimiter','\n'); % read all lines
header_long=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(1:end),'[\t]'); % split coloumns
Table=CellTable2StrTable(Table_); % convert Cell table to string table
GnNames=Table(:,1);
INDEX=zeros(size(DatabaseNames,1),1);
SPLEN=zeros(size(DatabaseNames,1),1);
Proteins=size(GnNames);

for i=1:1:Proteins
    found=find(strcmp(DatabaseNames,GnNames{i}));
    if(isempty(found)==0)
        INDEX(found)=1;
        SPLEN(found)=str2double(Table{i,5});
    end
end
INDEX=logical(INDEX);
end

function [Threshold]=DefineScoreThres(Scores,groups,groups_cyto,groups_peripheral,title_,figNo)
%% Distribution of prediction scores
legend_names={'Secreted' 'F1' 'Cytoplasmic'};
No_SECR=sum(groups);
No_CYTO=sum(groups_cyto);
No_PERIPH=length(groups)-No_SECR-No_CYTO;
MAX_Scores=max(Scores);
%% Otsu threshold
Scores_=(Scores./MAX_Scores);
[Threshold_R]=Otsu(Scores_(or(groups,groups_peripheral)),No_SECR+No_CYTO,'h');
[Threshold_L]=Otsu(Scores_(or(groups_cyto,groups_peripheral)),No_SECR+No_CYTO,'h');
Threshold_L=(Threshold_L*MAX_Scores);
Threshold_R=(Threshold_R*MAX_Scores);
%% Relative frequency
step=max(Scores);
x=min(Scores):(max(Scores)-min(Scores))./10:max(Scores);
[freq_Scores_secr]=hist(Scores(groups),x);freq_Scores_secr=freq_Scores_secr./No_SECR;
[freq_Scores_cyto]=hist(Scores(groups_cyto),x);freq_Scores_cyto=freq_Scores_cyto./No_CYTO;
[freq_Scores_periph]=hist(Scores(groups_peripheral),x);freq_Scores_periph=freq_Scores_periph./No_PERIPH;

%% Figure titles and legends
[h1]=FigureLegends(x,[freq_Scores_secr' freq_Scores_periph' freq_Scores_cyto'],figNo,'Score','Percent',title_,legend_names,'b');
hold on;
y_=0:.05:max([freq_Scores_secr freq_Scores_cyto freq_Scores_periph]);
x_L=ones(length(y_),1)*Threshold_L;
x_R=ones(length(y_),1)*Threshold_R;
plot(x_L,y_,'r.-',x_R,y_,'ro-');
text(Threshold_L,y_(length(y_)),'A|F1','Color','k','FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','bottom');
text(Threshold_L,y_(length(y_)),num2str(Threshold_L),'Color','k','FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top');
text(Threshold_R,y_(length(y_)),'F1|SEC','Color','k','FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','bottom');
text(Threshold_R,y_(length(y_)),num2str(Threshold_R),'Color','k','FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top');
Threshold=[Threshold_L Threshold_R];
hold off;
%% Save figures
% saveas(h1,['Figures/Patches_',int2str(windows(1)),'-',int2str(windows(end)),'.bmp'],'bmp');

end