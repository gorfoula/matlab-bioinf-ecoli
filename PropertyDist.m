function []=PropertyDist(fig,PropSum,x_label,fig_title)

colors_=['m','b','y','k','c','o'];
files=max(PropSum(:,1));

mean_val=zeros(files,1);
f_handle=figure(fig);
for i=1:1:files
    figure(fig);hold on;
    data=PropSum(PropSum(:,1)==i,2);
    [map x]=hist(data);
    bar(x,(map./length(data))*100,colors_(i));
    mean_val(i)=mean(data);
end
legend('Peri (SP len)','Cyto','B-barrel','Location','Best');
xlabel(x_label);
ylabel('Percent');
title(fig_title);

y_axe=1:4:50;
for i=1:1:files
    x_axe=mean_val(i)*ones(length(y_axe),1);
    plot(x_axe,y_axe,'o','MarkerEdgeColor','k','MarkerFaceColor',colors_(i));
end
s1 = regexp(fig_title,'[()]');
scale=fig_title(s1(1)+1:s1(2)-1);
pos=fig_title(s1(end-1):s1(end));
saveas(f_handle,['Figures/',x_label,'_',scale,'_',pos,'.bmp'],'bmp');

end