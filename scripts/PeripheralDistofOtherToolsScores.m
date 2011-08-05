function []=PeripheralDistofOtherToolsScores()
close all;
[CATCH] = ReadTable('JoinFiles/k12/(DB)TableS3_K12_ProteomeV13(for)ProteomeV13_BasicProteome_SEQ(1)F.txt','\n');
[TABLE Address]=Douplicates('JoinFiles/Peripheral/LIT_ECOLI.txt',[0 1],CATCH);

CATCH_=CATCH;

% MaskConv
hydro_len=19;
window=1:hydro_len;
num_windows=5;
table=[];count=0;
Y=[];groups=[];Y2=[];groups2=[];Y3=[];groups3=[];Y4=[];groups4=[];
mid=(hydro_len-1)/2;
thres=[19 35];

ind=find(strcmp(TABLE(:,1),'B') | strcmp(TABLE(:,1),'F1') | strcmp(TABLE(:,1),'A'));  % 19-37 285-309 447-469 554-570

for i=ind'
    cur_prob_tab=CATCH_(Address==i,:);
    table=[];amphi_table=[];amphi_len=[];
    table_counts=zeros(size(cur_prob_tab,1),1);
    table_mean_profile=zeros(size(cur_prob_tab,1),1);
    for j=1:1:size(cur_prob_tab,1)
        if(length(cur_prob_tab{j,3})>=hydro_len)
            [ProfileDist trash amphi]=ScaleProfile('K',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'ai',0);
            [ProfileK PhoScaleK]=ScaleProfile('K',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'g',0);
            ProfileC=ScaleProfile('C',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'g',0);
            Profile1=ScaleProfile('WW1',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'g',0);
            Profile2=ScaleProfile('WW2',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'g',0);
            Profiled=ScaleProfile('WWd',cur_prob_tab{j,2},cur_prob_tab{j,3},hydro_len,'g',0);
                       
            Profile1=-Profile1;
            Profile2=-Profile2;            
            Profiled=-Profiled;
            
            Profile_tm=Profiled;
            Profile=.2*ProfileK+.75*Profile_tm+.05*ProfileC;
%             Profile=0.4*(.5*ProfileK+.5*ProfileC)+0.6*Profile_tm;
            norm_pos=mean(amphi,2)./length(ProfileK);
            ncterm=norm_pos<.10 | norm_pos>.90;
            amphi_table=[amphi_table;sum(ncterm)];
            amphi_len=[amphi_len;amphi(ncterm,2)-amphi(ncterm,1)];
            
      
%             if(exist('cf','var')==1);close(cf);end
%             cf=figure;th=[];
%             [cf]=PlotAll(Profile1,Profile,ProfileK,Profiled,ProfileC,ProfileDist,amphi,cf,cur_prob_tab(j,:));
            
            count_patches=0;
%             pos=1;
%             while (pos<=length(Profile))
%                 
%                 n1=pos;
%                 n2=pos;
%                 patch_len=0;
%                 cur_max=Profile(pos);
%                 while(pos+1<=length(Profile) & Profile(pos+1)>=0 & Profile(pos+1)-cur_max>-.01)
%                     patch_len=patch_len+1;
%                     pos=pos+1;
%                     n2=pos;
%                     cur_max=Profile(pos);
%                 end
%                 cur_min=cur_max;
%                 while(pos+1<=length(Profile) & Profile(pos+1)>=0 & cur_min-Profile(pos+1)>-0.01)
%                     patch_len=patch_len+1;
%                     pos=pos+1;
%                     n2=pos;
%                     cur_min=Profile(pos);
%                 end
%                 if(patch_len==0)
%                     pos=pos+1;
%                 else
% %                     if(sum(Profile_tm(n1:n2)>-0.02)./(n2-n1+1)>.2)
% 
%                             count_patches=count_patches+floor(patch_len/hydro_len);
%                             table=[table;[n1 n2 patch_len] trapz(n1:n2,Profile(n1:n2))];
%                                                       
% %                             figure(cf);hold on;
% %                             plot(n1:n2,ones(1,n2-n1+1),'r','LineWidth',13);
% %                             text((n1+n2)/2,1,num2str(patch_len),'VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
% % %                             text((n1+n2)/2,0,num2str(sum(Profile_tm(n1:n2)>0)./(n2-n1+1)),'VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
% %                             text(n1,0.95,num2str(n1),'Color','r','VerticalAlignment','top');
% %                             text(n2,0.95,num2str(n2),'Color','r','VerticalAlignment','top');
% %                             hold off;
% %                             
% %                     end
%                 end
%                 table_counts(j)=count_patches;
%             end
        end
    end
%     Y=[Y;table(:,end)];
%     groups=[groups;ones(length(table),1)*count];
% 
%     Y2=[Y2;table(:,end-1)];
%     groups2=[groups2;ones(length(table),1)*count];
    
    Y3=[Y3;amphi_table];
    groups3=[groups3;ones(length(amphi_table),1)*count];
    
    Y4=[Y4;amphi_len];
    groups4=[groups4;ones(length(amphi_len),1)*count];

    count=count+1;
end
categoty=TABLE(ind,1);
x=0:1:max(Y3);
FigureLegends(x,Y3,groups3,[{''} {'# of a-helical amphiphilic  TMs'} {'Percent over all a-helical amphiphilic TMs'}],categoty,'p',{'-','.';':','o';'-','.'},0,[1 2 1 1]);


categoty=TABLE(ind,1);
x=min(Y4):2:max(Y4);
FigureLegends(x,Y4,groups4,[{''} {'Len of a-helical amphiphilic  TMs'} {'Percent over all a-helical amphiphilic TMs'}],categoty,'p',{'-','.';':','o';'-','.'},0,[1 2 1 2]);


% categoty=TABLE(ind,1);
% x=min(table(:,4)):(max(table(:,4))-min(table(:,4)))/50:max(table(:,4));
% FigureLegends(x,Y,groups,[{''} {'AUC of a-helical TMs'} {'Percent over all a-helical TMs'}],categoty,'p',{'-','.';':','o';'-','.'},0);
% % FigureLegends(1:max(Y),Y,groups,[{''} {'Hydrophobic Island Lenght'} {'Percent over all islands'}],categoty,'p',{'-','.';':','o';'-','.'},0);
% % FigureLegends((0:0.1:1)*max(Y2),Y2,groups2,[{''} {['Num of hydrophobic islands (',num2str(thres(1)),'-',num2str(thres(2)),')']} {'Percent of proteins'}],categoty,'p',{'-','.';':','o';'-','.'},0);
% FigureLegends((0:0.1:1)*max(Y2),Y2,groups2,[{''} {'a-helical TMs Lenght'} {'Percent over all a-helical TMs'}],categoty,'p',{'-','.';':','o';'-','.'},0);
% 
% x=0:.01:1;
% x=x*max(CellTable2Double(CATCH_(:,6)));
% 
% Cytoplasmic=CellTable2Double(CATCH_(Address==2,6));
% F1=CellTable2Double(CATCH_(Address==3,6));
% B=CellTable2Double(CATCH_(Address==5,6));
% Y=[Cytoplasmic;F1;B];
% groups=[zeros(length(Cytoplasmic),1);ones(length(F1),1);ones(length(B),1)*2];
% FigureLegends(x,Y,groups,[{''} {'Exp number of AAs in TMHs in n-term'} {'Percent of proteins'}],[TABLE(3,1) TABLE(2,1) TABLE(5,1)],'p',{'-','o';'-','.';':','x'},0);
% 
% Cytoplasmic=CellTable2Double(CATCH_(Address==2,end));
% F1=CellTable2Double(CATCH_(Address==3,end));
% B=CellTable2Double(CATCH_(Address==5,end));
% Y=[Cytoplasmic;F1;B];
% groups=[zeros(length(Cytoplasmic),1);ones(length(F1),1);ones(length(B),1)*2];
% FigureLegends(x,Y,groups,[{''} {'Probability of N-term to be in cyto side'} {'Percent of proteins'}],[TABLE(3,1) TABLE(2,1) TABLE(5,1)],'p',{'-','o';'-','.';':','x'},0);


end


function [cf]=PlotAll(Profile1,Profile,ProfileK,Profiled,ProfileC,ProfileDist,amphi,cf,cur_prob_tab)

x=1:length(Profile);
zero_axis=zeros(length(x),1);
figure(cf);subplot(2,2,1:2);
h=plot(x,ProfileK,'r','MarkerSize',1);
h=plot(x,Profile,'y',x,ProfileK,':c',x,Profiled,':k',x,ProfileC,'r',x,zero_axis,'g',x,ProfileDist,'b:o','MarkerSize',1);
%%%% print predicted amphiphilic  %%%
for cnt=1:size(amphi,1)
    if (exist('th','var') & isempty(th)==0);EraseText(th);end                
    figure(cf);subplot(2,2,1:2);hold on;
    plot(amphi(cnt,1):amphi(cnt,2),ones(1,amphi(cnt,2)-amphi(cnt,1)+1)*0.2,'r','LineWidth',10);
    text((amphi(cnt,1)+amphi(cnt,2))/2,0.2,num2str(amphi(cnt,2)-amphi(cnt,1)+1),'VerticalAlignment','middle','HorizontalAlignment','center');
    hold off;
    figure(cf);subplot(2,2,3);[th]=DrawHelix(cur_prob_tab{3}(amphi(cnt,1):amphi(cnt,2)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(cf);subplot(2,2,1:2);title(cur_prob_tab{2}(1:10));
figure(cf);subplot(2,2,1:2);
ylim([-1 1]);
% hold on; plot(1:length(ProfileK),ones(1:length(ProfileK),1)*.2,':r'); hold off;
legend('Kyte+(DGoct-DGifw)+Connet','Kyte','DGoct-DGifw','Connet(a-helix)');
            
cell_str=Str2CellArray(cur_prob_tab{3});
index=strcmp(cell_str,'H');

set(gca,'XTick',find(index));
set(gca,'XTickLabel',MergeColumns([cell_str(index) Double2CellTable(find(index),1)],''));
set(gca,'FontSize',8);
            
if(sum(index)>0)
    hold on;h = stem(find(index),Profile1(index),'Color','b');
    set(h,'Marker','o','MarkerFaceColor','k','LineStyle','none');hold off;
end
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% print experimental amphipath helixes  %%
AmphiNames={'GLPD_ECOLI' 'DNAA_ECOLI' 'MIND_ECOLI' 'DHNA_ECOLI' 'PTGA_ECOLI' 'PBPB_ECOLI'};
n1_=[357 355 259 390 1 284];n2_=[374 370 270 406 18 300];

loci=find(CellTable2Double(regexp(cur_prob_tab{2},AmphiNames)));
if(isempty(loci)==0)
    figure(cf);subplot(2,2,1:2);hold on;
    plot(n1_(loci):n2_(loci),zeros(1,n2_(loci)-n1_(loci)+1),'c','LineWidth',13);
    text((n1_(loci)+n2_(loci))/2,0,num2str(n2_(loci)-n1_(loci)+1),'VerticalAlignment','middle','HorizontalAlignment','center','FontWeight','bold');
    text(n1_(loci),-0.15,num2str(n1_(loci)),'Color','c','VerticalAlignment','top');
    text(n2_(loci),-0.15,num2str(n2_(loci)),'Color','c','VerticalAlignment','top');
    figure(cf);subplot(2,2,4);DrawHelix(cur_prob_tab{3}(n1_(loci):n2_(loci)));
    hold off;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end