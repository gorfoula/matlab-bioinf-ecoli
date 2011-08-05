function []=AminoDistrFig(Values,Names,save_,outFile)

% Values:       a matrix with sequence/feature representation
% Names:        Names of Representation (Categories/amino acids)
% save_:        1 save figures
%               0 not save figure
% outFile:      out file path of figures

[Row,Col,Zax]=size(Values);

rand('twister', sum(100*clock));
R=rand(54,1);
G=0:.01:1;
B=rand(length(R),1);
CATEGORIES=length(Names);

for i=1:1:Col
    f_handle=figure(2);
    x=1:1:CATEGORIES;
    [map,x_tr]=hist(Values(:,i),x);
    percent=(map./Row).*100;
    bar(x,percent);  %%%  precentage of all proteins having this property/amino acid
    title(['Position: ',num2str(i)]);
    xlabel('Amino acid');
    ylabel('Presentage over all (%)');
    for catg=1:1:CATEGORIES;
        text(catg,percent(catg),Names{catg},'Color',[R(catg) G(catg) B(catg)],'FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
    end
    grid;
    pause;
    if(save_==1)
        saveas(f_handle,[outFile,int2str(i),'.jpg'],'jpg');
    end 
end

close all;


%     text(2,map(2),'A','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(3,map(3),'V','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(4,map(4),'L','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(5,map(5),'I','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(6,map(6),'P','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(7,map(7),'F','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(8,map(8),'W','Color','y','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(9,map(9),'G','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(10,map(10),'S','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(11,map(11),'C','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(12,map(12),'N','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(13,map(13),'Q','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(14,map(14),'Y','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(15,map(15),'T','Color','k','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(16,map(16),'K','Color','r','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(17,map(17),'R','Color','r','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(18,map(18),'H','Color','r','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(19,map(19),'D','Color','b','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
%     text(20,map(20),'E','Color','b','FontSize',10,'FontWeight','demi','HorizontalAlignment','center','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k');
