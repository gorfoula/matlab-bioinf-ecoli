function []=SFGraph(inputFILE,Table,FEATURES)

%%%%%%%%%%%%%%%%%%  READ_FILES   %%%%%%%%%%%%%%%%%%%%%
%[aa ID NC PPos SecrMean CytoMean AllMean] = textread(inputFILE','%f %d %d %d %s %f %f %f',-1,'delimiter','\t');

[FEATURES_color]=AssignFeatureColor(length(FEATURES));
% 
% FEATURES_color(1,:)=[0 1 0];    %V,L,I,M,F green
% FEATURES_color(2,:)=[1 1 0];
% FEATURES_color(3,:)=[0 0 1];
% FEATURES_color(4,:)=[0 0 1];
% FEATURES_color(5,:)=[1 0 0];
% FEATURES_color(6,:)=[1 0 0];
% FEATURES_color(7,:)=([150 144 144]./255);
% FEATURES_color(8,:)=[0 0 0];    %V,L,I,M,F green
% FEATURES_color(9,:)=([255 128 0]./256);
% FEATURES_color(10,:)=([200 128 0]./256);
% FEATURES_color(11,:)=[1 0 1];

Sorted=sortrows(Table,[3 -4 -5]);

aa=Sorted(:,1);
NC=Sorted(:,2);
PPos=Sorted(:,3);
SecrMean=Sorted(:,4);
CytoMean=Sorted(:,5);
AllMean=Sorted(:,6);

OLC=FEATURES(NC);
OLC_pos1=OLC(PPos<=2);
max_numletters=0;
for i=1:1:length(OLC_pos1)
    max_numletters=max(max_numletters,length(OLC_pos1{i}));
end

fact=10;
aa=(aa*fact);
% aa=log2(100./aa);
% aa=(aa*100)./max(aa);

% aa(SecrMean<CytoMean)=(-1)*aa(SecrMean<CytoMean);
mean_dif=SecrMean;
mean_dif(aa<=0)=CytoMean(aa<=0);%1-SecrMean(aa<=0);
mean_dif(aa>0)=SecrMean(aa>0);

% SECRETED=2./(-log2(mean_dif));
SECRETED=mean_dif*6;

max_pos=max(PPos);
min_pos=min(PPos);
[freq,x]=hist(PPos,[min_pos:1:max_pos]);

[Wsum]=SumWeigthsPerPosition(PPos,aa,min_pos,max_pos);

f_handle=figure(1);
bar([min_pos:1:max_pos],Wsum,'w','BarWidth',1,'LineWidth',1);

% Max_fig=max([max(aa)+max(SECRETED) max(Wsum)]);
% Min_fig=min([min(aa)-min(SECRETED) min(Wsum)]);
Max_fig=max(Wsum);
Min_fig=min(aa)-min(SECRETED);

ylim([Min_fig-(fact) Max_fig+(fact)]);
xlim([min(PPos)-4 100+4]);

% ylim([Min_fig-Max_letter_size Max_fig+Max_letter_size]);
% xlim([min(PPos)-Max_letter_size max(PPos)+4+Max_hor_letter_size]);


% s1 = regexp(inputFILE, '[\/]');
% s2 = regexp(inputFILE, '[_.]');
% s3 = regexp(inputFILE, '[.]');
% name=inputFILE(s1(end)+1:s2(1)-1);
% savename=inputFILE(s1(end)+1:s3(1)-1);
% 
% grid on;
% xlabel('Position on Sequence');
% ylabel(['Sum of feature weights (x',int2str(fact),')']);
% title(name);


for i=1:1:length(SecrMean)
    i
PPos(i)
aa(i)
OLC(i)
FEATURES_color(NC(i),:)
abs(SECRETED(i))
    if(aa(i)<0)    
        text(PPos(i),aa(i),OLC(i),'Color',FEATURES_color(NC(i),:),'FontUnits','centimeters','FontSize',abs(SECRETED(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
    elseif(aa(i)>0)
        text(PPos(i),aa(i),OLC(i),'Color',FEATURES_color(NC(i),:),'FontUnits','centimeters','FontSize',abs(SECRETED(i)),'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end


% saveas(f_handle,['Figures/',savename,'.bmp'],'bmp');


end