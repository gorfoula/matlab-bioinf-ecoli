
% [CATCH]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/k12/PredictionsOfTools/EcoliK12_tmhmm_test.txt',[2 1],[1 1],[5 32]);
% CATCH_=CATCH(2:end,:);
% [TABLE Address]=Douplicates('JoinFiles/k12/PredictionsOfTools/EcoliK12_tmhmm.txt',[0 1],CATCH_);
% FastaSubSelection('JoinFiles/k12/TableS3_K12_ProteomeV13_BasicProteome.txt','JoinFiles/k12/ECOLI_K12_4407.fasta','mn',[1 1]);

clear all;
close all;

% [HEAD SEQ]=FastaRead('JoinFiles/k12/ECOLI_K12_4407.TableS3_K12_ProteomeV13_BasicProteome.fasta');
% FileWriteTable('JoinFiles/k12/ProteomeV13_BasicProteome_SEQ.txt',[HEAD' SEQ'],[],'w');

% [CATCH]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/k12/ProteomeV13_BasicProteome_SEQ.txt',[1 1],[1 0],[5]);
[CATCH] = ReadTable('JoinFiles/k12/(DB)TableS3_K12_ProteomeV13(for)ProteomeV13_BasicProteome_SEQ(1)F.txt','\n');
[TABLE Address]=Douplicates('JoinFiles/Peripheral/LIT_ECOLI.txt',[1 1],CATCH);

CATCH_=CATCH(2:end,:);

% MaskConv
hydro_len=19;
window=1:hydro_len;
num_windows=5;
table=[];count=0;
Y=[];groups=[];
Y2=[];groups2=[];
mid=(hydro_len-1)/2;
thres=[6 50];

ind=find(strcmp(TABLE(:,1),'F1') | strcmp(TABLE(:,1),'B') | strcmp(TABLE(:,1),'A'));

for i=ind'
    cur_prob_tab=CATCH_(Address==i,:);
    table=[];
    table_counts=zeros(size(cur_prob_tab,1),1);
    table_mean_profile=zeros(size(cur_prob_tab,1),1);
    for j=1:1:size(cur_prob_tab,1)
        if(length(cur_prob_tab{j,3})>=hydro_len)
            
            ProfileK=ScaleProfile('K',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,0);
            Profile1=ScaleProfile('WW1',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,0);
            Profile2=ScaleProfile('WW2',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,0);
            Profiled=ScaleProfile('WWd',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,0);
            
            x=1:length(Profile1);
            zero_axis=zeros(length(x),1);
            h=plot(x,ProfileK,'c',x,Profile1,'y',x,Profile2,'r',x,Profiled,'k',x,zero_axis,'g');
            legend('Kyte','DGwoct','DGwif','DGwoct-DGwif');
            
            cell_str=str2CellArray(cur_prob_tab{j,3});
            index=strcmp(cell_str,'H');
% 
%             set(gca,'XTick',find(index));
%             set(gca,'XTickLabel',MergeColumns([cell_str(index) Double2CellTable(find(index),1)],''));
%             set(gca,'FontSize',8);
            
            if(sum(index)>0)
                hold on;h = stem(find(index),Profile1(index),'Color','b');
                set(h,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold off;
            end
            
            Profile1=-Profile1;
            Profile2=-Profile2;            
            Profile=-Profiled;
            
            Profile=Profile2;
           
            count_patches=0;
            pos=1;
            while (pos<=length(Profile))
                
                n1=pos;
                n2=pos;
                patch_len=0;
                while(pos<=length(Profile) & Profile(pos)>0)
                    patch_len=patch_len+1;
                    pos=pos+1;
                    n2=pos;
                end
                if(patch_len==0)
                    pos=pos+1;
                else
                    if(patch_len<=thres(2) && patch_len>=thres(1))
                        table=[table;[n1 n2 patch_len]];
                        count_patches=count_patches+1;
                    end
                end
                table_counts(j)=count_patches;
            end
        end
    end
    Y=[Y;table(:,end)];
    groups=[groups;ones(length(table),1)*count];

    Y2=[Y2;table_mean_profile];
    groups2=[groups2;ones(length(table_counts),1)*count];

    count=count+1;
end
categoty=TABLE(ind,1);
FigureLegends(1:max(Y),Y,groups,[{''} {'Hydrophobic Island Lenght'} {'Percent over all islands'}],categoty,'p',{'-','.';':','o';'-','.'},0);
FigureLegends((0:0.1:1)*max(Y2),Y2,groups2,[{''} {['Num of hydrophobic islands (',num2str(thres(1)),'-',num2str(thres(2)),')']} {'Percent of proteins'}],categoty,'p',{'-','.';':','o';'-','.'},0);

x=0:.01:1;
x=x*max(CellTable2Double(CATCH_(:,6)));

Cytoplasmic=CellTable2Double(CATCH_(Address==2,6));
F1=CellTable2Double(CATCH_(Address==3,6));
B=CellTable2Double(CATCH_(Address==5,6));
Y=[Cytoplasmic;F1;B];
groups=[zeros(length(Cytoplasmic),1);ones(length(F1),1);ones(length(B),1)*2];
FigureLegends(x,Y,groups,[{''} {'Exp number of AAs in TMHs in n-term'} {'Percent of proteins'}],[TABLE(3,1) TABLE(2,1) TABLE(5,1)],'p',{'-','o';'-','.';':','x'},0);

Cytoplasmic=CellTable2Double(CATCH_(Address==2,end));
F1=CellTable2Double(CATCH_(Address==3,end));
B=CellTable2Double(CATCH_(Address==5,end));
Y=[Cytoplasmic;F1;B];
groups=[zeros(length(Cytoplasmic),1);ones(length(F1),1);ones(length(B),1)*2];
FigureLegends(x,Y,groups,[{''} {'Probability of N-term to be in cyto side'} {'Percent of proteins'}],[TABLE(3,1) TABLE(2,1) TABLE(5,1)],'p',{'-','o';'-','.';':','x'},0);