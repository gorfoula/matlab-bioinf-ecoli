function []=SignalPeptideMatureFusions(IFNames,sel)
%% INPUT:
%%  IFNames:    Cell table of file names
close all;
outDir='Gems/Curated/';
ModelNames={'CytoVSsecr';'CytoVSsecr_mat_FIN'};
%% Read Database File
[text_] = textread(IFNames{1},'%s',-1,'delimiter','\n'); % read all lines
header_long=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(1:end),'[\t]'); % split coloumns
Table_Secr=CellTable2StrTable(Table_); % convert Cell table to string table
%% 
[text_] = textread(IFNames{2},'%s',-1,'delimiter','\n'); % read all lines
header_long=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_(1:end),'[\t]'); % split coloumns
Table_Cyto=CellTable2StrTable(Table_); % convert Cell table to string table

Colors=RainbowColor(80);
% col_index=1:floor(length(Colors)./categories):length(Colors);
% Colors=Colors(col_index,:);

SEQ_cyto=Table_Cyto(:,6);
SEQ=Table_Secr(:,6);
SPLen=CellTable2Double(Table_Secr(:,5));
SP_Names=Table_Secr(:,1);

MATS=size(SEQ_cyto,1);
SIGNALS=size(SEQ,1);

%%              SIGNALS
%%       ___ ___ ___ ___ ___
%%      |___|___|___|___|___|
%%      |___|___|___|___|___|       MATURES
%%      |___|___|___|___|___|
%%
%%      ROW=mod(Pos,SIGNALS)
%%      COL=floor(Pos/SIGNALS)+(row>0)

SPSeqComb=cell(SIGNALS,SIGNALS);
MatSeqComb=cell(MATS,SIGNALS);
SPLenComb=zeros(SIGNALS,SIGNALS);
MatLenComb=zeros(MATS,SIGNALS);
for peptide=1:1:MATS
    SEL_MATURE=SEQ_cyto{peptide}(2:end);     % Methionine excluded
    for i=1:1:SIGNALS
        if(peptide==1750)
            display(peptide);
        end
        MatSeqComb{peptide,i}=[SEQ{i}(1:SPLen(i)),SEL_MATURE];
        MatLenComb(peptide,i)=length(MatSeqComb{peptide,i});
        if(peptide<=SIGNALS)
            SPSeqComb{peptide,i}=[SEQ{i}(1:SPLen(i)),SEQ{peptide}(SPLen(peptide)+1:end)];
            SPLenComb(peptide,i)=length(SPSeqComb{peptide,i});
        end
    end
end

SPLen=SPLen-1;
mask=zeros(MATS,SIGNALS);mask(:,1)=1;  %% create a mask in order to replicate SPLen
SPLen2D=conv2(SPLen',mask);

if(sel==1)  %% Load Scores if they exist
    [SCORES_sp(:,1) SCORES_sp(:,2) SCORES_sp(:,3)] = textread('Data/Shuffles/SPxSP_Scores.txt','%f %f %f',-1,'delimiter','\t');
    [SCORES_mat(:,1) SCORES_mat(:,2) SCORES_mat(:,3)] = textread('Data/Shuffles/SPxMAT_Scores.txt','%f %f %f',-1,'delimiter','\t');
    [INDEX_Len_sp] = textread('Data/Shuffles/SPxSP_INDEX.txt','%d',-1,'delimiter','\t');INDEX_Len_sp=logical(INDEX_Len_sp);
    [INDEX_Len_mat] = textread('Data/Shuffles/SPxMAT_INDEX.txt','%d',-1,'delimiter','\t');INDEX_Len_mat=logical(INDEX_Len_mat);
else
    [SCORES_sp INDEX_Len_sp]=ComboPredictor(ModelNames,reshape(SPSeqComb,SIGNALS*SIGNALS,1),reshape(SPLen2D(1:SIGNALS,1:SIGNALS),SIGNALS*SIGNALS,1),reshape(SPLenComb,SIGNALS*SIGNALS,1),outDir,0);
    [SCORES_mat INDEX_Len_mat]=ComboPredictor(ModelNames,reshape(MatSeqComb,MATS*SIGNALS,1),reshape(SPLen2D(1:MATS,1:SIGNALS),MATS*SIGNALS,1),reshape(MatLenComb,MATS*SIGNALS,1),outDir,0);

    dlmwrite('Data/Shuffles/SPxSP_Scores.txt', SCORES_sp, '\t')
    dlmwrite('Data/Shuffles/SPxMAT_Scores.txt',SCORES_mat, '\t')
    dlmwrite('Data/Shuffles/SPxSP_INDEX.txt',INDEX_Len_sp, '\t')
    dlmwrite('Data/Shuffles/SPxMAT_INDEX.txt',INDEX_Len_mat, '\t')
end

%% Score Dist
% [Scores_norm]=NormalizeScores(SCORES_sp(:,3));

[BvsW_index]=CalculateANDplot(SCORES_sp(:,1),SCORES_mat(:,1),INDEX_Len_sp,INDEX_Len_mat,SP_Names,SIGNALS,MATS,1);
PlotFeatureDiff(ModelNames,SPSeqComb,SPLenComb,SPLen,SP_Names,outDir,BvsW_index,SIGNALS);
% CalculateANDplot(SCORES_sp(:,2),SCORES_mat(:,2),INDEX_Len_sp,INDEX_Len_mat,SP_Names,SIGNALS,MATS,4);
% CalculateANDplot(SCORES_sp(:,3),SCORES_mat(:,3),INDEX_Len_sp,INDEX_Len_mat,SP_Names,SIGNALS,MATS,5);
end

function []=PlotFeatureDiff(ModelNames,SPSeqComb,SPLenComb,SPLen,SP_Names,outDir,BvsW_index,SIGNALS)

for j=1:SIGNALS
    [SCORES_sp INDEX_Len_sp Code_best W]=ComboPredictor(ModelNames,SPSeqComb(BvsW_index(j,1),j),SPLen(j),SPLenComb(BvsW_index(j,1),j),outDir,0);
    [SCORES_sp_worst INDEX_Len_sp Code_worst]=ComboPredictor(ModelNames,SPSeqComb(BvsW_index(j,2),j),SPLen(j),SPLenComb(BvsW_index(j,2),j),outDir,0);
    
    if(SCORES_sp(1)>0 && SCORES_sp_worst(1)<0)
        index_n=strcmp({'PPB' 'GLNH' 'OMPA' 'MACA' 'NRFB'},SP_Names{j}(1:end-6));
        if(sum(index_n)>0)
            VALUES=20;
            FEATURES={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};
        
            temp1=zeros(VALUES*100,1);temp1(Code_best(:,1))=1;
            temp2=zeros(VALUES*100,1);temp2(Code_worst(:,1))=1;
            res=find(and(logical(temp1)==1,logical(temp2)==0));
    
            NumericCode=mod(res,VALUES);
            PosOnPepseq=floor(res/VALUES)+(NumericCode>0);
            NumericCode(NumericCode==0)=VALUES;
            OneLetterCode=FEATURES(NumericCode);
            DrawFeatures(23,W(and(logical(temp1)==1,logical(temp2)==0)),PosOnPepseq,OneLetterCode,NumericCode,VALUES);
    
            res=find(and(logical(temp1)==0,logical(temp2)==1));
    
            NumericCode=mod(res,VALUES);
            PosOnPepseq=floor(res/VALUES)+(NumericCode>0);
            NumericCode(NumericCode==0)=VALUES;
            OneLetterCode=FEATURES(NumericCode);
    
            hold on;
            display(SP_Names{j});
            DrawFeatures(23,W(and(logical(temp1)==0,logical(temp2)==1)),PosOnPepseq,OneLetterCode,NumericCode,VALUES,'red');   
            figure(23);hold on;
            title([SP_Names{j}(1:end-6),' BEST/WORST: ',SP_Names{BvsW_index(j,1)}(1:end-6),' / ',SP_Names{BvsW_index(j,2)}(1:end-6),' ',num2str(SCORES_sp(1)),' / ',num2str(SCORES_sp_worst(1))]);
            hold off;
            pause;
        end
    end

end

end

function [BestVSworst_index]=CalculateANDplot(Scores_sp,Scores_mat,INDEX_Len_sp,INDEX_Len_mat,SP_Names,SIGNALS,MATS,fig)
x=(-1:.1:1)*max(abs([Scores_sp;Scores_mat]));x=[-2 1.5 3]+1;
Scores_sp_WT=Scores_sp(and(INDEX_Len_sp,reshape(eye(SIGNALS),SIGNALS*SIGNALS,1)>0));
[dist]=hist(Scores_sp,x);dist=(dist*100)./(SIGNALS*SIGNALS);
[dist_WT]=hist(Scores_sp_WT,x);dist_WT=(dist_WT*100)./(SIGNALS);
[dist_mat]=hist(Scores_mat,x);dist_mat=(dist_mat*100)./(SIGNALS);

[Pos]=find(Scores_sp==3);Pos=Pos(1);ROW=mod(Pos,SIGNALS);COL=floor(Pos/SIGNALS)+(ROW>0);
SPTable=reshape(Scores_sp,SIGNALS,SIGNALS);

indx=find(strcmp('PPB_ECOLI',SP_Names));
x=(-1:0.1:1)*max(abs(SPTable(indx,:)));
[freq]=hist(SPTable(indx,:),x);
bar(x,(freq*100)./sum(freq));
hold on; stem(SPTable(indx,indx),max((freq*100)./sum(freq)),'.-r');hold off;

xlabel('Secretion Score');
ylabel('Percent over all combinations');

[best_score best_combo]=max(SPTable(indx,:));
x=(-1:0.1:1)*max(abs(SPTable(best_combo,:)));
[freq]=hist(SPTable(best_combo,:),x);
figure;
bar(x,(freq*100)./sum(freq));
hold on; stem(best_score,max((freq*100)./sum(freq)),'.-r');hold off;

xlabel('Secretion Score');
ylabel('Percent over all combinations');
title(SP_Names(best_combo));

[dist_BEST_FUSION]=hist(SPTable(:,COL),x);dist_BEST_FUSION=(dist_BEST_FUSION*100)./(SIGNALS);
[dist_PPB_MATURE]=hist(SPTable(find(strcmp(SP_Names,'PPB_ECOLI'),1),:),x);dist_PPB_MATURE=(dist_PPB_MATURE*100)./(SIGNALS);
[h5]=FigureLegends([1 2 3],dist_BEST_FUSION',fig+5,'Scores (aas)','Percent','Signal Peptide Shuffle',{},'b');
[h6]=FigureLegends([1 2 3],dist_PPB_MATURE',fig+6,'Scores (aas)','Percent','Signal Peptide Shuffle',{},'b');
[h7]=FigureLegends([1 2 3],dist',fig+7,'Scores (aas)','Percent','Signal Peptide Shuffle',{},'b');

display(dist_BEST_FUSION);
display(dist_PPB_MATURE);
display(dist);

[BEST_dist WORST_dist]=BESTandWorstCOMBO(Scores_sp,INDEX_Len_sp,SIGNALS);
[BEST_dist_mat WORST_dist_mat]=BESTandWorstCOMBO(Scores_mat,INDEX_Len_mat,SIGNALS);

BEST=sortrows([BEST_dist' (1:1:SIGNALS)'],[-1]);
WORST=sortrows([WORST_dist' (1:1:SIGNALS)'],[-1]);
BEST_mat=sortrows([BEST_dist_mat' (1:1:SIGNALS)'],[-1]);
WORST_mat=sortrows([WORST_dist_mat' (1:1:SIGNALS)'],[-1]);

legend_names={'SP and Mature shuffling';'WT';'Cyto mature and SP Combos'};

% Scores_sp_TEMP=reshape(Scores_sp,SIGNALS,SIGNALS);
% index_n=find(strcmp(SP_Names,{'PPB_ECOLI'}));
% [ppb_dist x]=hist(Scores_sp_TEMP(:,index_n),[-1 1 3]);
% [h1]=FigureLegends(x',(ppb_dist*100)'./SIGNALS,70,'Scores (aas)','Percent','PPB_ECOLI',[],'b');
% figure(70);hold on;
% y=0:5:100;
% x=ones(1,size(y))*Scores_sp_TEMP(index_n,index_n);
% plot(x,y,'r.-');


%% Figure titles and legends
[h1]=FigureLegends(x,[(dist'./(SIGNALS*SIGNALS)) (dist_WT'./SIGNALS) (dist_mat'./(SIGNALS*MATS))],fig,'Scores (aas)','Percent','Signal Peptide Shuffle',legend_names,'b');
[h2]=FigureLegends((1:1:SIGNALS),[BEST_dist' WORST_dist'],fig+1,'Signal Peptide','Percent','BEST and WORST Signal Peptide',{'BEST' 'WORST'},'pm',{':' '*';'-' 'o'});hold on;
GREY=Grayscale(2);
stem(BEST(1:7,2),BEST(1:7,1),'Color',GREY(1,:),'Marker','*');stem(WORST(1:7,2),WORST(1:7,1),'Color',GREY(2,:),'Marker','o');
for j=1:7   %for all best signal peptides
    text(BEST(j,2),BEST(j,1),SP_Names{BEST(j,2)}(1:end-6),'Rotation',25,'Color',GREY(1,:),'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    text(WORST(j,2),WORST(j,1),SP_Names{WORST(j,2)}(1:end-6),'Rotation',25,'Color',GREY(2,:),'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');  
end
hold off;
SCORES=reshape(Scores_sp,SIGNALS,SIGNALS);
BestVSworst_index=zeros(SIGNALS,4);
% [INDEX SPLEN]=MatchNames(SP_Names,'Dataset/AllExperimental.txt')
for j=1:SIGNALS   %for all best signal peptides
%     if(strcmp(SP_Names{j},'PPBKL_ECOLI'))
%     figure(fig+2);
%     stem(1:SIGNALS,SCORES(j,:));
%     title(SP_Names{j});
%     cur_scores=sortrows([SCORES(j,:)' (1:SIGNALS)'],[-1]); 
%     for i=1:7 %i=find(INDEX)'
% %         text(i,SCORES(j,i),SP_Names{i}(1:end-6),'Color','k','FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
%         text(cur_scores(i,2),cur_scores(i,1),SP_Names{cur_scores(i,2)}(1:end-6),'Color','k','FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
%         text(cur_scores(end-i+1,2),cur_scores(end-i+1,1),SP_Names{cur_scores(end-i+1,2)}(1:end-6),'Color','r','FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
%     end
        [val index]=max(SCORES(:,j));
        [val_min index_min]=min(SCORES(:,j));
        BestVSworst_index(j,:)=[index index_min val val_min];
%         SP_Names{index}
%         SP_Names{index_min}
        FileWriteTable('Data/Best_Worst_combos.txt',[SP_Names(index) num2str(index) SP_Names(index_min) num2str(index_min)],'Best Name\tScore\tWorst Name\tScore','a');
%         pause;
%     end
end
hold off;
[h3]=FigureLegends((1:1:SIGNALS),[BEST_dist_mat' WORST_dist_mat'],fig+3,'Signal Peptide','Percent','Cyto and SP combos',{'BEST' 'WORST'},'pm',{':' '*';'-' 'o'});hold on;
stem(BEST_mat(1:7,2),BEST_mat(1:7,1),'Color',GREY(1,:),'Marker','*');stem(WORST_mat(1:7,2),WORST_mat(1:7,1),'Color',GREY(2,:),'Marker','o');
for j=1:7   %for all best signal peptides
    text(BEST_mat(j,2),BEST_mat(j,1),SP_Names{BEST_mat(j,2)}(1:end-6),'Rotation',25,'Color',GREY(1,:),'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    text(WORST_mat(j,2),WORST_mat(j,1),SP_Names{WORST_mat(j,2)}(1:end-6),'Rotation',25,'Color',GREY(2,:),'FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');  
end
hold off;
% %% Save figures
% saveas(h1,'Figures/SPShuffling_ScoreDist_.eps','eps');
% saveas(h2,'Figures/BEST_WORST_SPs_.eps','eps');
% saveas(h3,'Figures/Cyto_SP_Compos_.eps','eps');


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
    found=find(strcmpi(DatabaseNames,[GnNames{i},'_ECOLI']));
    if(isempty(found)==0)
        INDEX(found)=1;
        SPLEN(found)=str2double(Table{i,5});
    end
end
INDEX=logical(INDEX);
end

function [BEST WORST]=BESTandWorstCOMBO(Scores,INDEX_Len,SIGNALS)
%% Choose BEST and WORST Signal Peptide
VAR=std(Scores);  % variance of scores
MEAN=mean(Scores); % mean value of scores
% INDEX_WORST=find(and(INDEX_Len,Scores<0));
% INDEX_BEST=find(and(INDEX_Len,Scores>0));
INDEX_WORST=find(and(INDEX_Len,Scores<(-(2*VAR)+MEAN)));
INDEX_BEST=find(and(INDEX_Len,Scores>(2*VAR+MEAN)));


BEST_SP=mod(INDEX_BEST,SIGNALS);BEST_SP(BEST_SP==0)=SIGNALS;
[BEST_dist]=hist(BEST_SP,1:1:SIGNALS);
BEST=BEST_dist*100./SIGNALS;

WORST_SP=mod(INDEX_WORST,SIGNALS);WORST_SP(WORST_SP==0)=SIGNALS;
[WORST_dist]=hist(WORST_SP,1:1:SIGNALS);
WORST=WORST_dist*100./SIGNALS;
end