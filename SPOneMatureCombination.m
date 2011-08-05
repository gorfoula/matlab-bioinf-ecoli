function []=SPOneMatureCombination(IFile,CytoFile,SelFeatFILE,WeigthFile,len,sel)

close all;
Red=rand(54,1);Green=0:.01:1;Blue=rand(length(Red),1); %%%  random colors for plots
%%% Read Secreted and Cytoplasmic Files
[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[GnNames_cyto, Descr_cyto, TotLen_cyto, MW_cyto, SPLen_cyto, SEQ_cyto] = textread(CytoFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read cyto dataset file
%%% Read Selected Features File
[aa, FeatureID] = textread(SelFeatFILE,'%d %d',-1,'delimiter','\t');
FeatureID=FeatureID-1;
[w] = textread(WeigthFile,'%f',-1,'delimiter','\t');  %%% read weights

VALUES=20;  %%% division factor / number of attribute IDs
if(exist('len','var')==0)
    NumericCode=mod(max(FeatureID),VALUES);
    PosOnPepseq=floor(max(FeatureID)/VALUES)+(NumericCode>0);
    len=PosOnPepseq;
end

%%% SECRETED Only Peptides with sufficient amino acid length
MOLECULARW=MW(TotLen>=SPLen+len);
DESCRIPTION=Descr(TotLen>=SPLen+len);
Names=GnNames(TotLen>=SPLen+len);   %Gene Names
SELECTED=SEQ(TotLen>=SPLen+len);    %Select Petides with Seq long enough
SIGNALS=length(SELECTED);
TOTALLEN=TotLen(TotLen>=SPLen+len);
SPLen=SPLen(TotLen>=SPLen+len);     %Select corresponding Signal Peptide Lengths

%%% CYTOPLASMIC Only Peptides with sufficient amino acid length
sel_index=TotLen_cyto>=SPLen_cyto+len;
MOLECULARW_cyto=MW_cyto(sel_index);
DESCRIPTION_cyto=Descr_cyto(sel_index);
Names_cyto=GnNames_cyto(sel_index);   %Gene Names
SELECTED_cyto=SEQ_cyto(sel_index);    %Select Petides with Seq long enough
SIGNALS_cyto=length(SELECTED_cyto);
MATURES=length(SELECTED_cyto);
TOTALLEN_cyto=TotLen_cyto(sel_index);
SPLen_cyto=SPLen_cyto(sel_index);     %Select corresponding Signal Peptide Lengths

[allscores allindex WT_SP_score]=CalculateScores(FeatureID,w,Names,SELECTED,SELECTED,SPLen,SPLen,VALUES,len,sel);MATURES=SIGNALS;
% [allscores allindex WT_SP_score]=CalculateScores(FeatureID,w,Names,SELECTED,SELECTED_cyto,SPLen,SPLen_cyto,VALUES,len,sel);

% best_SP_index=reshape(best_SP_index,1,BEST*samples);
% best_SP_score=reshape(best_SP_score,1,BEST*samples);
allscores_row=reshape(allscores,1,MATURES*SIGNALS);
allindex_row=reshape(allindex,1,MATURES*SIGNALS);

[thres]=PlotTotalScoreHistogram(allscores_row,WT_SP_score,MATURES)
best_SP_score=allscores_row(allscores_row>=thres(1));
worst_SP_score=allscores_row(allscores_row<=thres(2));

best_SP_index=allindex_row(allscores_row>=thres(1));
worst_SP_index=allindex_row(allscores_row<=thres(2));

[Sorted_index]=PrintBestSPHist(SelFeatFILE,best_SP_index,Names,SIGNALS,'BestSP',3);
[Sorted__worst_index]=PrintBestSPHist(SelFeatFILE,worst_SP_index,Names,SIGNALS,'WorstSP',4);
WritebestSPGroupfiles(best_SP_index,Sorted_index,Names,DESCRIPTION,TOTALLEN,MOLECULARW,SPLen,SELECTED,SIGNALS,'best');
WritebestSPGroupfiles(worst_SP_index,Sorted__worst_index,Names,DESCRIPTION,TOTALLEN,MOLECULARW,SPLen,SELECTED,SIGNALS,'worst');
% PrintScoreCurvesEachGroup(allscores,Names,Sorted_index,best_SP_index,samples)

end

function [allscores allindex WT_SP_score]=CalculateScores(FeatureID,w,Names,SP_pool,MAT_pool,SPLen_SP_pool,SPLen_MAT_pool,VALUES,len,sel)
SIGNALS=length(SP_pool);
MATURES=length(MAT_pool);
WT_SP_score=zeros(MATURES,1);
% best_SP_index=zeros(samples,BEST);
% best_SP_score=zeros(samples,BEST);
allscores=zeros(MATURES,SIGNALS); % all combination scores
allindex=zeros(MATURES,SIGNALS);  % index of what?
index_exp=FindExperimental(Names);
index_exp=index_exp(index_exp>0);
phoAKL_score=0;
for peptide=1:1:MATURES
    cur_MATURE=MAT_pool{peptide}(SPLen_MAT_pool(peptide)+1:end);
    AAs3D=zeros(SIGNALS,VALUES,len);
    SelSeq=cell(SIGNALS,1);
    cur_Sequence=cell(SIGNALS,1);
    for i=1:1:SIGNALS    % transfrom to 0-1 coding mode  ALL SP Combinations with specific MATURE
        cur_Sequence{i}=[SP_pool{i}(2:SPLen_SP_pool(i)),cur_MATURE];  %% for on mature and all Signal Peptides
        [AAs3D(i,:,:) SelSeq(i)]=AARepresentation('bin',1,len,{cur_MATURE});  %% Methionine excluded
    end
    
    [CatgNames,Catg3D]=CatgRepr(AAs3D,'bin','Features/mature_AAProperties.txt');
    [CatgNames_9f,Catg3D_9f]=CatgRepr(AAs3D,'bin','Features/AAProperties_new.txt');
    
    AAs3D=reshape(AAs3D,SIGNALS,VALUES*len);
    Catg3D=reshape(Catg3D,SIGNALS,len*11);
    Catg3D_9f=reshape(Catg3D_9f,SIGNALS,len*9);
    
    if(sel==20119)
        AAs3D=[AAs3D Catg3D Catg3D_9f];
    end

    score=AAs3D(:,FeatureID)*w(1:end-1)+w(end);aa=1:1:SIGNALS;
    allscores(peptide,:)=score;
    allindex(peptide,:)=1:1:SIGNALS;

%     sorted_score=sortrows([aa' score],[-2]);
%     best_SP_index(peptide,:)=sorted_score(1:BEST,1);
%     best_SP_score(peptide,:)=sorted_score(1:BEST,2);
    if(MATURES==SIGNALS)
        WT_SP_score(peptide)=score(peptide);
    end
    
%     if(or(strcmpi(Names{peptide},'phoAKL'),strcmpi(Names{peptide},'phoA')))        %%% print phoA and phoAKL scores
%         f6=figure(6);hold on;
%         if(strcmpi(Names{peptide},'phoAKL'))
%             h = stem(score,'fill','--');
%             set(get(h,'BaseLine'),'LineStyle',':');
%             set(h,'MarkerFaceColor','c');
%             color='k';
%             phoAKL_score=score;
%         else
%             stem(score);color='r';
%             score_error=mean(sqrt((score-phoAKL_score).^2))
%             title(['Score PhoA VS phoAKL, MSQE : ',num2str(score_error)]);
%             xlabel('Signal Peptide');
%             ylabel('Score');
%             legend(legend('phoAKL','phoA'));
%             
%             sum(sign(score(index_exp)).*sign(phoAKL_score(index_exp))<0)
%             temp=(sign(score).*sign(phoAKL_score)<0).*(1:1:length(score))';
%             temp=temp(index_exp);
%             temp=temp(temp>0);
%             {Names{temp}}
%         end
%         
%         for sp=index_exp     % again for all best signal peptides
%             text(sp,score(sp),Names{sp},'Color',color,'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
%         end
%         hold off;
%         saveas(f6,['Figures/PhoA_PhoAKL_score.bmp'],'bmp');
% 
%     end
end

end

function [thres]=PlotTotalScoreHistogram(allscores_row,WT_SP_score,samples)

f2=figure(2);hold on;
[m s]=normfit(allscores_row);
thres=[2.5*s+m m-2.5*s];
[freq x]=hist(WT_SP_score,min(allscores_row):1:max(allscores_row));
[freq_all x]=hist(allscores_row,min(allscores_row):1:max(allscores_row));hold on;
bar(x,(freq_all*100)./length(allscores_row),'m');
bar(x,(freq*100)./samples,0.4,'g');
y_vertical=0:1:max((freq_all*100)./length(allscores_row));
x_down_thr=ones(length(y_vertical))*thres(1);
x_up_thr=ones(length(y_vertical))*thres(2);
plot(x_down_thr,y_vertical,x_up_thr,y_vertical);
title('Score hist');
legend('All SP Mature Combinations','Wild type SP Mature combination','Location','best');
hold off;
saveas(f2,['Figures/WT_SPs_score_hist.bmp'],'bmp');
end

function [Sorted_index]=PrintBestSPHist(SelFeatFILE,best_SP_index,Names,samples,fname,fig_count)
Red=rand(54,1);Green=0:.01:1;Blue=rand(length(Red),1);
f3=figure(fig_count);               % hist of all combination scores
[freq x]=hist(best_SP_index,1:1:samples);
freq=(freq*100)./length(best_SP_index);
bar(x,freq);hold on;
index=find(freq>0);
freq_=freq(freq>0);
Sorted_index=sortrows([index' freq_'],[-2]);

count=1;
for j=Sorted_index(:,1)'   %for all best signal peptides
    if(Sorted_index(count,2)>=1.2)    
        figure(fig_count);hold on;  % print name of best signal peptide  
        text(x(j),freq(j),Names{j},'Color',[Red(count) Green(count) Blue(count)],'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
        hold off;
    end
    count=count+1;
end
title(['SPs that give the',fname,'score for Secreted Proteins']);
ylabel('Percentage over all Proteins');
found=regexp(SelFeatFILE,'[/.]');
hold off;
saveas(f3,['Figures/',fname,'(',SelFeatFILE(found(end-1)+1:found(end)-1),').bmp'],'bmp');
end

function []=WritebestSPGroupfiles(SP_index,Sorted_index,Names,DESCRIPTION,TOTALLEN,MOLECULARW,SPLen,SELECTED,samples,type)
Red=rand(54,1);Green=0:.01:1;Blue=rand(length(Red),1);
[freq x]=hist(SP_index,1:1:samples);
count=1;
s_best_idex=fopen(['Data/All',type,'SPIndex.txt'],'w');
for j=Sorted_index(:,1)'   %for all best signal peptides
    group=find(SP_index==x(j));
%     if(j==157)
%         display('Im lost without you');
%     end
%     if(Sorted_index(count,2)>=1.2)
%         s=fopen(['Data/',Names{j},'_',type,'_SP_with.txt'],'w');  %% file save group of proteins with specific best SP
%         fprintf(s,'Names \t Description \t Total Len \t MW \t SP Len \t Sequence \t SP Sequence \t MATURE Sequence\n');
%         for ig=group   %for all protein with the same best signal peptide
%             if(ig==407)
%                 display('Im lost without you also.');
%             end
% 
%             combined=[SELECTED{j}(1:SPLen(j)),SELECTED{ig}(SPLen(ig)+1:end)];
%             fprintf(s,'%s \t %s \t %d \t %d \t %d \t %s \t %s \t %s \n',Names{ig},DESCRIPTION{ig},TOTALLEN(ig),MOLECULARW(ig),SPLen(ig),combined,SELECTED{ig}(1:SPLen(ig)),SELECTED{ig}(SPLen(ig)+1:end));
%         end
%         fclose(s);
%     end
    count=count+1;
    fprintf(s_best_idex,'%s \t %s \t %d \t %s \n',Names{j},DESCRIPTION{j},SPLen(j),SELECTED{j}(1:SPLen(j)));  %% file save best SP in an index
        
end
fprintf(s_best_idex,'%d \n',round(mean(SPLen(Sorted_index(:,1)))));  %% print mean value of SPLen
fclose all;

end

function []=PrintScoreCurvesEachGroup(allscores,Names,Sorted_index,best_SP_index,samples)
Red=rand(54,1);Green=0:.01:1;Blue=rand(length(Red),1);
[freq x]=hist(best_SP_index,1:1:samples);
for j=Sorted_index(:,1)'   %for all best signal peptides
    group=find(best_SP_index==x(j));
    color=1;
    for best=Sorted_index(:,1)'     % again for all best signal peptides
        f5=figure(5);hold on;       % plot score curves for each SP in the current group of proteins
        [freq_score x_score]=hist(allscores(group,best));
        degree=3;
        p = polyfit(x_score,freq_score,degree); % Degree 3 fit
        y_= polyval(p,x_score);
        plot(x_score,y_);
        text(mean(allscores(group,best)),max(y_),Names{best},'Color',[Red(color) Green(color) Blue(color)],'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
        title(['Group of proteins with Best SP:',Names{j}]);
        xlabel('Score');
        ylabel('Percent over all proteins');
        hold off;
        color=color+4;
    end
    saveas(f5,['Figures/ScoreCurves_',Names{j},'_group.bmp'],'bmp');close(5);    
end
end


function [index]=FindExperimental(Names)

[GnNames_exp, Descr_exp, TotLen_exp, MW_exp, SPLen_exp, SEQ_exp] = textread('Dataset/AllExperimental.txt','%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

num_exp=length(GnNames_exp);
index=zeros(1,num_exp);

for i=1:1:num_exp
    found=find(strcmpi(Names,GnNames_exp{i}));
    if(isempty(found)==0)
        index(i)=found;
    end
end

end