%%%   normpdf, normcdf: cumulative function, erf: area
close all;
x=-4:0.1:4;
Norm_dist=normpdf(x,0,1);
figure(1);hold on;
z=1.96;
% for pos=z:.01:x(end)
%     y_2=0:.0001:normpdf(pos,0,1);
%     a_2=ones(length(y_2))*pos;
%     plot(a_2,y_2,'c');
% end
h=stem(z:0.001:x(end),normpdf(z:0.001:x(end),0,1),'c');
set(h(1),'Marker','.');
h=stem(x(1):0.001:-z,normpdf(x(1):0.001:-z,0,1),'c');
set(h(1),'Marker','.')

text(z+1,normpdf(z,0,1),'a/2','Color','k','FontUnits','points','FontSize',20,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
text(-(z+1),normpdf(z,0,1),'a/2','Color','k','FontUnits','points','FontSize',20,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
text(z,0,'z_a_/_2','Color','k','FontUnits','points','FontSize',20,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
text(-z,0,'- z_a_/_2','Color','k','FontUnits','points','FontSize',20,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');

plot(x,Norm_dist);hold off;
