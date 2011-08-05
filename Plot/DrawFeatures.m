function []=DrawFeatures(fig,W,POS,LETTERS,NC,FTRS,edgcol)

letters=length(LETTERS);
[FEATURES_color]=AssignFeatureColor(FTRS);
figure(fig);
x_temp=-5:.1:max(POS);
plot(x_temp,zeros(1,length(x_temp)),'g');
xlim([-5 max(POS)]);
ylim([-5 5]);
if(exist('edgcol','var')==0)
    edgcol='none';
end
stem(POS(W~=0),W(W~=0),'k')
for i=1:1:letters
    if(W(i)<0)    
        text(POS(i),W(i),[LETTERS{i},'(',num2str(POS(i)),')'],'Rotation',45,'EdgeColor',edgcol,'Color',FEATURES_color(NC(i),:),'FontWeight','bold','HorizontalAlignment','Right','VerticalAlignment','top');
    elseif(W(i)>0)
        text(POS(i),W(i),[LETTERS{i},' (',num2str(POS(i)),')'],'Rotation',45,'EdgeColor',edgcol,'Color',FEATURES_color(NC(i),:),'FontWeight','bold','HorizontalAlignment','Left','VerticalAlignment','bottom');
    end
    hold on;
end
hold off;
end 