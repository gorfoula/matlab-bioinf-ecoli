function []=ComboPrediction(type)
% DBIFile,ModelNames,TestSet,sel,gems,wSelFiles,outDir
outDir='Gems/Curated/';
ModelNames={'CytoVSsecr';'CytoVSsecr_mat'};

TestSet_lit={'Data/LIT_MUT_secr.txt' 'Data/LIT_MUT_cyto.txt'};
TestSet={'Gems/Curated/TEST/AllSecr_test.txt' 'Gems/Curated/TEST/AllCyto_test.txt'};
TrainSet={'Gems/Curated/TRAIN/AllSecr.txt' 'Gems/Curated/TRAIN/AllCyto.txt'};
AllSet={'Gems/Curated/AllSecr.txt' 'Gems/Curated/AllCyto.txt'};
MYCTUSet={'Data/MYCTU_valSet_Leversen_mn.txt' 'Data/MYCTU_valSet_Leversen_NEG_mn.txt'};

DBCurated='Data/DataSetGO/DB_SeptV1.txt';
MutantDB='Data/LIT_MUT_mn_ToolPred.txt';
DBMYCTU='Data/MYCTU_valSet_Leversen_all_mn_ToolPred.txt';

[SCORES GROUPS]=Prediction(DBCurated,ModelNames,AllSet,0,1,0,outDir);

Peptides=size(GROUPS,1);
AUC_all=[];
switch (type)
    case 'acc'
        cp = classperf(double(GROUPS));
        classes=zeros(Peptides,1);
        classes(or(SCORES(:,1)>=0,and(SCORES(:,1)<0,SCORES(:,2)>=0)))=1;
        classperf(cp,classes,[1:1:Peptides]');
        display(['<Accurancy> COMBO MODEL: ',num2str(cp.CorrectRate)]);
    case 'auc'
        close all;
        %% best combo model
        x_perc=0:.01:1;
        for perc=x_perc
            allscores_=(perc*SCORES(:,1))+((1-perc)*SCORES(:,2));
            AUC_=AUC(GROUPS, NormalizeScores(allscores_), 1);
            AUC_all=[AUC_all;AUC_];
        end
        %% FIGURE  %%%%
        [h1]=FigureLegends(x_perc,AUC_all,11,'Weigth of Model 1 (including SP)','AUC','AUC of combined Model',[],'p');
        hold on;
        best=max(AUC_all);
        y_=AUC_all;
        x_max=ones(length(y_),1)*x_perc(AUC_all==best);
        plot(x_max,y_,'r.-');
        text(best,y_(end),['BEST AUC: ',num2str(best),' at ',num2str(x_max(1))],'Color','k','FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','bottom');
        hold off;
        %%
        display(['BEST AUC: ',num2str(best)]);
        display(['RULE AUC: ',num2str(AUC_all(x_perc==.45))]);
        %% Accuracy of best combo
        perc=x_max(1);
        allscores_=(perc*SCORES(:,1))+((1-perc)*SCORES(:,2));
        display(['<Accurancy> COMBO MODEL: ',num2str(Accuracy(GROUPS,allscores_))]);
        scores_combo=(0.45*SCORES(:,1))+((1-0.45)*SCORES(:,2));
        display(['<Accurancy> COMBO MODEL: ',num2str(Accuracy(GROUPS,scores_combo))]);
end
end