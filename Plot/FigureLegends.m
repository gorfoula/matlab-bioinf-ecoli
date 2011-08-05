function [Colors]=FigureLegends(x,Y,groups,Label,legend_,type,linetype,grey_,splot)
% close all;
%% if type='pm' or (bm) the mean on y-axis is ploted
%% if type='p.-' or (bm) the mean on y-axis is ploted
%% Colors
if(exist('splot','var')==0)      %if expr is not defined
    figure;
else
    figure(splot(1));
    if(length(splot)>1)
        subplot(splot(2),splot(3),splot(4));
    end
end
hold on;
n=max(groups)+1;
[y x mx mn]=PDF(Y,groups,x);
if(grey_==0) %%%%%%   COLORS   %%%%%%%%%%%%%%%%
    Colors=RainbowColor(100);
    col_index=1:floor(length(Colors)./n):length(Colors);
    Colors=Colors(col_index,:);
else         %%%%%%   GREY   %%%%%%%%%%%%%%%%
    GREY=Grayscale(n);
    Colors=GREY(1:end-1,:);
end
colormap(Colors);
% if(n>2)
% % %     [p,h] = ttest2(truey(:,1),truey(:,2));  % Mann Whitney test
% %     [p,h] = ranksum(y(:,1),y(:,2));  % Mann Whitney test
% %     title_=[title_,' Mann Whitney: ',num2str(p),' ',num2str(h)];
%     [p1,h1,s1] = anova1(truey,groups,'off');     display(p1);
%     [t,m] = multcompare(s1); pause;
%     [p,h,s] = kruskalwallis(truey,groups,'off');     display(p);
%     [t,m] = multcompare(s);pause;    
%     title_=[title_,' ANOVA p-value: ',num2str(p1),' ','Kruskal Wallis p-value: ',num2str(p)];    
% end
%%%%%%%%%%%%%%%  PLOT DISTRIBUTIONS   %%%%%%%%%%%%%%%%%%%%%%
switch (type(1))
    case 'p'
        for i = 1:n
            if(isempty(linetype)==0)
%                 p = polyfit(x',y(:,i),2);
%                 f = polyval(p,x);
                if(strcmp(linetype{i,2},''))
%                     plot(x,f,'Color',Colors(i,:),'LineWidth',2,'linestyle',linetype{i,1});hold on;
                    plot(x,y(:,i),'Color',Colors(i,:),'LineWidth',2,'linestyle',linetype{i,1});hold on;
                else
%                     plot(x,y(:,i),'Color','k','LineWidth',4,'linestyle',linetype{i,1});hold on;
%                     plot(x,f,'Color',Colors(i,:),'LineWidth',1,'linestyle',linetype{i,1},'Marker',linetype{i,2},'MarkerSize',4,'MarkerFaceColor',Colors(i,:));hold on;
                    plot(x,y(:,i),'Color',Colors(i,:),'LineWidth',1.2,'linestyle',linetype{i,1},'Marker',linetype{i,2},'MarkerSize',3,'MarkerFaceColor',Colors(i,:));hold on;
                end
            else
                plot(x,y(:,i),'Color',Colors(i,:),'LineWidth',2);hold on;
            end
        end
    case 'b'
        h=bar(x,y,1,'grouped');
%         h=bar3(x,y,0.35,'detached');
        set(h,'LineWidth',1);
        for i = 1:n
            ch = get(h(i),'Children');
            fvd = get(ch,'Faces');
            fvcd = get(ch,'FaceVertexCData');
            fvcd(:) = i;
            set(ch,'FaceVertexCData',fvcd)
        end
    otherwise
        display('Wrong type of plot, p/b  fr plot/bar');
end
%%%%%%%%%%%%%%%%%%%%  PLOT MEAN VALUES  %%%%%%%%%%%%%
for rv=1:1:n %% for all random variables
    y_axis=(0:.05:1)*mx(rv);
    x_axis=ones(length(y_axis),1);
    if(strcmp(linetype{i,2},''))
        plot(x_axis*mn(rv),y_axis,'Color',Colors(rv,:),'LineWidth',1,'linestyle',linetype{rv,1});
    else
        plot(x_axis*mn(rv),y_axis,'Color',Colors(rv,:),'LineWidth',0.8,'linestyle',linetype{rv,1},'MarkerFaceColor',Colors(rv,:),'MarkerSize',3,'Marker',linetype{rv,2});
    end
    %%% grey [153 153 153]./256
    text(mn(rv),mx(rv),num2str(mn(rv)),'BackgroundColor',[1 1 1],'HorizontalAlignment','center','VerticalAlignment','baseline');
end
title(Label{1});
xlabel(Label{2});
ylabel(Label{3});
legend(legend_);
hold off;
end