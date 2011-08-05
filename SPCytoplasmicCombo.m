function []=SPCytoplasmicCombo(SecrFile,CytoFile,SelFeatFILE,WeigthFile,len)

Red=rand(54,1);
Green=0:.01:1;
Blue=rand(length(Red),1);
% close all;
[GnNames_cyto, Descr_cyto, TotLen_cyto, MW_cyto, SPLen_cyto, SEQ_cyto] = textread(CytoFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(SecrFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[aa, FeatureID] = textread(SelFeatFILE,'%d %d',-1,'delimiter','\t');
FeatureID=FeatureID-1;
VALUES=20;

COLORS={'r', 'g', 'b', 'c', 'm', 'y', 'k', 'w'};

if(exist('len','var')==0)
    NumericCode=mod(max(FeatureID),VALUES);
    PosOnPepseq=floor(max(FeatureID)/VALUES)+(NumericCode>0);
    len=PosOnPepseq;
end

SAMPLES=length(SEQ_cyto);
SIGNALS=length(SEQ);

best_SP_index=zeros(SAMPLES,1);
best_SP_score=zeros(SAMPLES,1);
allscores=ones(SAMPLES,SIGNALS)*(-100);
allscores_secr=ones(SIGNALS,SIGNALS)*(-100);



for peptide=1:1:SAMPLES
    SEL_MATURE=SEQ_cyto{peptide}(2:end);     % Methionine excluded
    AAs3D=zeros(SIGNALS,VALUES,len);
    SelSeq=cell(SIGNALS,1);
    AAs3D_secr=zeros(SIGNALS,VALUES,len);
    for i=1:1:SIGNALS
%         i
        if((SPLen(i)+TotLen_cyto(peptide)-2)>=len)
            cur_Sequence=[SEQ{i}(2:SPLen(i)),SEL_MATURE];
            [AAs3D(i,:,:) SelSeq(i)]=AARepresentation('bin',1,len,{cur_Sequence});  
        end
        if(peptide<=SIGNALS && (TotLen(peptide)-SPLen(peptide)>=len))
            cur_Sequence=[SEQ{i}(2:SPLen(i)),SEQ{peptide}(2:end)];
            [AAs3D_secr(i,:,:) SelSeq(i)]=AARepresentation('bin',1,len,{cur_Sequence});  
        end
    end

    AAs3D=reshape(AAs3D,SIGNALS,VALUES*len);
    
    [w] = textread(WeigthFile,'%f',-1,'delimiter','\t');
    score=AAs3D(:,FeatureID)*w(1:end-1)+w(end);
    allscores(peptide,:)=score;
    
    if(peptide<=SIGNALS)
        AAs3D_secr=reshape(AAs3D_secr,SIGNALS,VALUES*len);
    
        [w] = textread(WeigthFile,'%f',-1,'delimiter','\t');
        score_secr=AAs3D_secr(:,FeatureID)*w(1:end-1)+w(end);
        allscores_secr(peptide,:)=score_secr;
    end
    
    max_score=max(score);max_score_index=find(score==max_score);
    best_SP_index(peptide)=max_score_index(1);
    best_SP_score(peptide)=max_score;

end


allscores_row=reshape(allscores,1,SAMPLES*SIGNALS);  %%% all score of cyto and SP combinations in a row
allscores_secr_row=reshape(allscores_secr,1,SIGNALS*SIGNALS);   %%% all score of Secreted and SP combinations in a row
allscores_secr_row=allscores_secr_row(allscores_secr_row>-100);

[thres]=PlotTotalScoreHistogram(allscores_row,WT_SP_score,samples)
best_SP_score=allscores_row(allscores_row>=thres(1));
worst_SP_score=allscores_row(allscores_row<=thres(2));


f2=figure(2);
x_cyto=min(allscores_secr_row):1:max(allscores_secr_row);
[freq_all x]=hist(allscores_row,x_cyto);hold on;perc_cyto=(freq_all*100)./length(allscores_row);
bar(x_cyto,perc_cyto,'r');
[freq_secr x_secr]=hist(allscores_secr_row,x_cyto);hold on;perc_secr=(freq_secr*100)./length(allscores_secr_row);
bar(x_secr,perc_secr,0.4,'g');
    mean_cyto=ones(1,max(perc_secr) )*mean(allscores_row);
    y_cyto=1:1:max(perc_secr);
    plot(mean_cyto,y_cyto,'o','MarkerEdgeColor','k','MarkerFaceColor','g');
    mean_secr=ones(1,max(perc_secr) )*mean(allscores_secr_row);
    y_secr=1:1:max(perc_secr);
    plot(mean_secr,y_secr,'o','MarkerEdgeColor','k','MarkerFaceColor','r');
title('Score hist');
legend('All SP Cyto Combinations','All SP Secreted Mature Combinations','Location','best');
hold off;
f3=figure(3);
[freq x]=hist(best_SP_index,1:1:SIGNALS);
freq=(freq*100)./SAMPLES;
bar(x,freq);hold on;
index=find(freq>0);
freq_=freq(freq>0);
Sorted_index=sortrows([index' freq_'],[-2]);
count=1;
for j=Sorted_index(:,1)'   %for all best signal peptides
    group=find(best_SP_index==x(j));
    
    color=1;
    
    s=fopen(['Data/',GnNames{j},'_best_SP_with.txt'],'w');
    fprintf(s,'GnNames \t Description \t Total Len \t MW \t SP Len \t Sequence \t SP Sequence \t MATURE Sequence\n');
    for ig=group'   %for all protein with the same best signal peptide
        fprintf(s,'%s \t %s \t %d \t %d \t %d \t %s \t %s \t %s \n',GnNames_cyto{ig},Descr_cyto{ig},TotLen_cyto(ig),MW_cyto(ig),SPLen_cyto(ig),SEQ_cyto{ig});
    end
    fclose all;
    
    figure(3);hold on;
    text(x(j),freq(j),GnNames{j},'Color',[Red(count) Green(count) Blue(count)],'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    hold off;
    count=count+1;
    
end
title('SPs that give the best score for Secreted Proteins');
ylabel('Percentage over all Proteins');
hold off;




saveas(f2,['Figures/ScoreSPMAtureCombo_cytoVSsecr.bmp'],'bmp');
saveas(f3,['Figures/BestSPsHist_cyto.bmp'],'bmp');


end