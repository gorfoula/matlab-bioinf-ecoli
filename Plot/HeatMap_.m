function []=HeatMap(Atr,bins)
Categories=size(Atr,2);

MAX_x=-Inf;MAX_y=-Inf;
MIN_x=Inf;MIN_y=Inf;
for i=1:Categories
    MAX_x=max( [ MAX_x max(Atr{i}(:,1)) ] );
    MIN_x=min( [ MIN_x min(Atr{i}(:,1)) ] );

    MAX_y=max( [ MAX_y max(Atr{i}(:,2)) ] );
    MIN_y=min( [ MIN_y min(Atr{i}(:,2)) ] );
end

x_=MIN_x:(MAX_x-MIN_x)/(bins-1):MAX_x;
y_=MIN_y:( (MAX_y-MIN_y) / (length(x_)-1)):MAX_y;
temp=hot;
hot_inv=temp(end:-1:1,:);

for i=1:Categories
    Atr1=Atr{i}(:,1);
    Atr2=Atr{i}(:,2);
    [n x]= hist3([Atr1 Atr2],'Ctrs',{x_;y_});
    h=pcolor(x{1,2},x{1,1},n);
    set(h, 'zdata', ones(size(n))*(i-1)*10);
    colormap(hot_inv) % heat map 
    % grid on;
    view(3);
    hold on;
end
colorbar;
hold off;
end