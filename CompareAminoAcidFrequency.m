function []=CompareAminoAcidFrequency(DatasetFiles,len,step)

% close all;
[FileNames] = textread(DatasetFiles,'%s',-1,'delimiter','\t');
AMINOS={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};
if(exist('len','var')==0)
    len=25;
end

    [occur_cyto distr_cyto]=AminoAcidFrequency(FileNames{1},len,step);
    [occur_secr distr_secr]=AminoAcidFrequency(FileNames{2},len,step);

M=max((max(occur_cyto)),(max(occur_secr)));
m=min((min(occur_cyto)),(min(occur_secr)));

CYTO=length(occur_cyto(:,1));
SECR=length(occur_secr(:,1));

diff_=zeros(20,ceil(len/step));

s=fopen(['Data/AAFreq_ttest_CytoVSsecr_',num2str(len),'.txt'],'w');
fprintf(s,'AA \t Same \t p-value\n');
f3=figure(3);
for j=1:1:20
    close all;
%     f1=figure(1);
%     [freq x]=hist(occur_cyto(:,j),m(j):1:M(j));dist=(freq./CYTO)*100;
%     bar(x,dist,'r');
%     hold on;
%     [freq_2 x_2]=hist(occur_secr(:,j),m(j):1:M(j));dist_2=(freq_2./SECR)*100;
%     [h p]=ttest2(freq_2,freq);h=not(h);
%     bar(x_2,dist_2,0.4);
%     xlabel(['Occurences on the first ',int2str(len),' aas of Mature']);
%     ylabel('Precentage over all Peptides');
%     title(['Histogram of number of <',AMINOS{j},'> occurences on the first ',int2str(len),' aas of Mature (ttest:',int2str(h),', p-value:',num2str(p),')']);
%     legend('Cytoplasmic','Secreted','Location','best');
%     saveas(f1,['Figures/Hist_occurences_(',AMINOS{j},')','_(',int2str(len),')','.bmp'],'bmp');
%     fprintf(s,'%s \t %d \t %1.3f \n',AMINOS{j},h,p);
    degree=7;
    extra=sum(distr_secr(j,:)==0);   %%% zero occurances
    zero_index=find(distr_secr(j,:)==0);
    distr_secr(j,zero_index)=1;

%     f2=figure(2);
    x=step:step:len;
    perc_cyto=(distr_cyto(j,:)./CYTO);
    perc_secr=(distr_secr(j,:)./(SECR+extra));
%     bar(x,perc_cyto,'r');hold on;
%     bar(x,perc_secr,0.4);
%     title(['Histogram of number of <',AMINOS{j},'> occurences per position over all proteins']);
%     legend('Cytoplasmic','Secreted','Location','best');
%     p = polyfit(x,perc_cyto,4); % Degree 4 fit
%     y_cyto=polyval(p,x);
%     p = polyfit(x,perc_secr,4); % Degree 4 fit
%     y_secr=polyval(p,x);
%     plot(x,y_cyto,'r-o',x,y_secr,'b-o');hold off;
%     saveas(f2,['Figures/PerPos_occurences_(',AMINOS{j},').bmp'],'bmp');
    
%     for degree=1:1:30
        
        ratio=log2(perc_secr./perc_cyto);
        [p S mu]= polyfit(x,ratio,degree); % fit polynomyal
        x__=(x-mu(1))./mu(2);
        y_ratio=polyval(p,x__);
        diff_(j,:)=y_ratio;
%         figure(3);
%         bar(x,ratio,'y');hold on;
%         plot(x,y_ratio,[COLORS{degree},'-']);
%         title( [ 'degree',num2str(degree),'error',num2str(sum((y_ratio-ratio).^2)) ] );
%         pause;
%     end

%     pause;
    

 end

index=[4 5 16 17 18 19 20];
count=1;
y_limits=[min(min(diff_)) max(max(diff_))];
for i=index
    figure(3);
    subplot(length(index),1,count);
    diff=diff_(i,:);
    bar(x,diff,'r');
    xlim([1 len]);
    ylim(y_limits);
    text(-1,max(diff),AMINOS{i},'Color','k','FontUnits','points','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    count=count+1;
end
    saveas(f3,['Figures/PerPos_Cyto_secr_diff_.bmp'],'bmp');
index=[9 10 11 12 13 14 15];
count=1;
for i=index
    f4=figure(4);
    subplot(length(index),1,count);
    diff=diff_(i,:);
    bar(x,diff,'r');
    xlim([1 len]);
    ylim(y_limits);
    text(-1,max(diff),AMINOS{i},'Color','k','FontUnits','points','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    count=count+1;
end
    saveas(f4,['Figures/PerPos_Cyto_secr_diff_2.bmp'],'bmp');
index=[1 2 3 6 7 8];
count=1;
for i=index
    f5=figure(5);
    subplot(length(index),1,count);
    diff=diff_(i,:);
    bar(x,diff,'r');
    xlim([1 len]);
    ylim(y_limits);
    text(-1,max(diff),AMINOS{i},'Color','k','FontUnits','points','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    count=count+1;
end
    saveas(f5,['Figures/PerPos_Cyto_secr_diff_3.bmp'],'bmp');

fclose all;
end