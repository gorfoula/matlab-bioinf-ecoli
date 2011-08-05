function []=PrintValues(fig,counts,PropSum,x_label,fig_title)

mean_val=mean(PropSum);
y_axe=1:4:50;
x_axe=mean_val*ones(length(y_axe),1);
plot(x_axe,y_axe,'o','MarkerEdgeColor','k','MarkerFaceColor',colors(counts));

end