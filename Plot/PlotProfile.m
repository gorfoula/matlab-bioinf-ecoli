function []=PlotProfile(fig,Name,SQ,Profile)

len=length(SQ);step=len/2;
subps=floor(len/step)+double(mod(len,step)>0);
pos=1;MAX=abs(max(Profile));
for s=1:1:subps
    figure(fig);subplot(subps,1,s);
    bar(pos:1:step+pos-1,Profile(pos:pos+step-1),'y');
    xlim([pos-1 step+pos]);
    ylim([-MAX MAX]);
    for i=pos:1:step+pos-1
        text(i,0,SQ(i),'Color','k','FontName','Helvetica Narrow','FontUnits','centimeters','FontSize',0.3,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
    end
    pos=pos+step;
end
if(mod(len,step)~=0)
    figure(fig);subplot(subps,1,s+1);
    bar(pos:len,Profile(pos:len),'y');
    xlim([pos-1 step+pos]);
    ylim([-MAX MAX]);
    for i=pos:1:mod(len,step)+pos-1
        text(i,0,SQ(i),'Color','k','FontName','Helvetica Narrow','FontUnits','centimeters','FontSize',0.3,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
    end
end
figure(fig);subplot(subps,1,1);
title(Name);
% pause;

end