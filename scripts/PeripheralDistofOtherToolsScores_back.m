
% file='JoinFiles/k12/TableS3_K12_ProteomeV13.txt';
% 
% tools_indes=[];
% [namefile dir]=IsolateFileName({file});
% outfile=[dir{1},namefile{1},'_douplicates.txt'];
% IDS_intact=INTABLE(index(1)+1:end,index(2));
% IDS=INTABLE(index(1)+1:end,index(2));

clear all;
close all;

[CATCH]=FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/k12/PredictionsOfTools/EcoliK12_tmhmm_test.txt',[2 1],[1 1],[5 32]);
CATCH_=CATCH(2:end,:);
[TABLE Address]=Douplicates('JoinFiles/k12/PredictionsOfTools/EcoliK12_tmhmm.txt',[0 1],CATCH_);
% MaskConv
hydro_len=19;
window=1:hydro_len;
num_windows=5;
table=[];count=0;
Y=[];groups=[];
Y2=[];groups2=[];
mid=(hydro_len-1)/2;
thres=[2 6];
for i=[2 3 5]
    cur_prob_tab=CATCH_(Address==i,:);
    table=[];
    table_counts=zeros(size(cur_prob_tab,1),1);
    for j=1:1:size(cur_prob_tab,1)
        if(length(cur_prob_tab{j,2})>=hydro_len)
            display(cur_prob_tab{j,3});
%             [Profile]=ScaleProfile('K',cur_prob_tab{j,3},cur_prob_tab{j,2},hydro_len,0);
            [Profile1]=ScaleProfile('WW1',cur_prob_tab{j,3},cur_prob_tab{j,2},hydro_len,0);
            [Profile2]=ScaleProfile('WW2',cur_prob_tab{j,3},cur_prob_tab{j,2},hydro_len,0);
            [Profiled]=ScaleProfile('WWd',cur_prob_tab{j,3},cur_prob_tab{j,2},hydro_len,0);
            
            x=1:length(Profile1);
%             plot(x,Profile1,'.-y',x,Profile2,'.r',x,Profiled,'.k');
            
            Profile1=-Profile1;
            Profile2=-Profile2;            
            Profile=-Profiled;

            hold on;plot(x,Profile1,'y',x,Profile2,'r',x,Profile,'k');hold off;
            Profile=Profile2;
            
%             if(length(cur_prob_tab{j,2})>=100)
%                 [te]=ScaleProfile('WW2',cur_prob_tab{j,3},cur_prob_tab{j,2}(1:100),hydro_len,1);
%             end
%             [d1] = gradient(double(Profile),1:1:length(Profile));
%             hold on;subplot(2,1,1);plot(1:1:length(Profile),d1,'k');hold off;
%             [d2] = gradient(d1,1:1:length(Profile));
%             hold on; subplot(2,1,1);plot(1:1:length(Profile),d2,'r');hold off;
%             hold on; subplot(2,1,1);h=stem(1:1:length(Profile),(abs(d1)<=.015 & Profile>0 & d2<0)*.1);set(h,'Color','k','LineWidth',2);hold off;
            
%             for wl=4:hydro_len
                wl=mid;
                end_=wl;
                count_patches=0;
                while(end_<=length(Profile))
                    cur_window=end_-wl+1:end_;
                    [m n]=max(Profile(cur_window));
                    n=(cur_window(1)-window(1))+n;
                    if(m>0)
                        if(n-mid>0 && n+mid<=length(Profile))
                            cond=Profile(n-mid:n+mid)>0;
%                     cond=Profile(n-mid:n+mid)<0 & Profile2(n-mid:n+mid)>Profile1(n-mid:n+mid);
                        end
                        if(n-mid<=0)
                            cond=Profile(1:n+mid)>0;
%                     cond=Profile(1:n+mid)<0 & Profile2(1:n+mid)>Profile1(1:n+mid);
                        end
                        if(n+mid>length(Profile))
                            cond=Profile(n-mid:end)>0;
%                     cond=Profile(n-mid:end)<0 & Profile2(n-mid:end)>Profile1(n-mid:end);
                        end
                        patch_len=sum(cond);
                        n1=find(cond,1,'first')+(cur_window(1)-window(1));
                        n2=find(cond,1,'last')+(cur_window(1)-window(1));
                        end_=n+wl+((wl-1)/2);
                        table=[table;[n1 n2 n patch_len]];
                        [n1 n2 n patch_len]

                        if(patch_len<=thres(2) && patch_len>=thres(1))
                            count_patches=count_patches+1;
                        end
                    else
                        end_=end_+wl;
                    end
                end
                table_counts(j)=count_patches;
%             end
        end
    end
    Y=[Y;table(:,2)];
    groups=[groups;ones(length(table),1)*count];

    Y2=[Y2;table_counts];
    groups2=[groups2;ones(length(table_counts),1)*count];

    count=count+1;
end
categoty=TABLE([2 3 5],1);
FigureLegends(1:max(Y),Y,groups,[{''} {'Hydrophobic Island Lenght'} {'Percent over all islands'}],categoty,'p',{'-','.';':','o';'-','.'},0);
FigureLegends(1:max(Y2),Y2,groups2,[{''} {['Num of hydrophobic islands (',num2str(thres(1)),'-',num2str(thres(2)),')']} {'Percent of proteins'}],categoty,'p',{'-','.';':','o';'-','.'},0);

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