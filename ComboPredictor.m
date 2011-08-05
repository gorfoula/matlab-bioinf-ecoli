function [SCORES INDEX_Len INDEX W]=ComboPredictor(Models,SEQ,SPLen,TotLen,outDir,sel)
%% Read Selected Features File
[aa_sp, FeatureID_sp] = textread([outDir,Models{1},'.txt'],'%d %d',-1,'delimiter','\t');
FeatureID_sp=FeatureID_sp-1;  %%%% we don't want first column which is for category
[aa, FeatureID_mat] = textread([outDir,Models{2},'.txt'],'%d %d',-1,'delimiter','\t');
FeatureID_mat=FeatureID_mat-1;  %%%% we don't want first column which is for category
%% Read weight file
[w_sp] = textread([outDir,Models{1},'.txtw.txt'],'%f',-1,'delimiter','\t');  %%% read weights
[w_mat] = textread([outDir,Models{2},'.txtw.txt'],'%f',-1,'delimiter','\t');  %%% read weights
VALUES=20;  %%% division factor / number of attribute IDs
len=100;

INDEX_Len=TotLen>=len;
allscores_sp=zeros(size(SEQ,1),1);
allscores_mat=zeros(size(SEQ,1),1);
[allscores_sp(INDEX_Len) Indx_sp W]=ScoreCalculation(SEQ(INDEX_Len),[],w_sp,FeatureID_sp,VALUES,len,sel);   %%%% Score calculation
[allscores_mat(INDEX_Len)]=ScoreCalculation(SEQ(INDEX_Len),SPLen(INDEX_Len),w_mat,FeatureID_mat,VALUES,len,sel);   %%%% Score calculation

scores_combo=(0.45*allscores_sp)+((1-0.45)*allscores_mat);

SCORES=[allscores_sp allscores_mat scores_combo];

INDEX=[Indx_sp'];

end