x=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 24 25 28 30 32 35 36 40 48 50 60 70 100 120 140];
y=[1 1 1 1 1 1 1 1 1 1 1.02 1.03 1.04 1.06 1.09 1.1 1.12 1.13 1.13 1.15 1.14 1.13 1.1 1.14 1.17 1.21 1.23 1.23 1.25 1.23 1.26 1.24 1.25 1.26 1.28 1.29 1.3 1.29 1.3 1.3];

x_54=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1 2 3 4 5 6 7 8 10 12 14 16 20 24 28 32 36 40];
y_54=[1 1 1 1 1 1 1.01 1.01 1.02 1.035 1.085 1.125  1.175 1.2 1.22 1.235 1.25 1.3 1.3 1.34 1.37 1.4 1.4 1.43 1.4 1.4 1.4];

x_51=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 24 28 32 40 48 56 64 72 80];
y_51=[1 1 1 1 1 1 1 1 1 1  1.01 1.02 1.02  1.0435 1.05 1.055 1.06 1.07  1.07 1.075 1.09 1.105  1.11 1.12 1.125 1.15  1.16 1.16 1.18 1.19 1.2 1.2 1.2 1.2]

y=y(x>=1);x=x(x>=1);
y_51=y_51(x_51>=1);x_51=x_51(x_51>=1);
y_54=y_54(x_54>=1);x_54=x_54(x_54>=1);

degree=4;

p = polyfit(x,y,degree) % Degree 3 fit
y_= polyval(p,x);

p = polyfit(x_51,y_51,degree) % Degree 3 fit
y_51_= polyval(p,x_51);

p = polyfit(x_54,y_54,degree) % Degree 3 fit
y_54_= polyval(p,x_54);



% figure(1);
% for i=1:1:4 
% p = polyfit(x_54,y_54,i) % Degree 3 fit
% y__= polyval(p,x_54);
% plot(x_54,y_54,'k',x_54,y__,'c');
% title(['degree: ',int2str(i)]);
% pause;
% end

order=2;  %%%% of derivative (dif )
len=0.5:1:2;
close all;

f2=figure(2);
subplot(2,2,1);
plot(x,y_);
deriv=diff(y_,order);
deriv1=diff(y_,1);
x_56_=x(1:end-2);
x_56_=x_56_(and(deriv1(1:end-1)>0.0099,abs(deriv)<=0.0099));
th1=ones(length(len))*min(x_56_);
th2=ones(length(len))*max(x_56_);
hold on;
title(['GBP Fluo 5.6, ',int2str(min(x_56_)),'-',int2str(max(x_56_))]);
stem(x(1:end-order),deriv);
find(deriv==0)

subplot(2,2,2);
deriv_54=diff(y_54_,order);
deriv1=diff(y_54_,1);
x_54_=x_54(1:end-2);
x_54_=x_54_(and(deriv1(1:end-1)>0.0099,abs(deriv_54)<=0.0099));
stem(x_54(1:end-order),deriv_54,'rs-');
th1_54=ones(length(len))*min(x_54_);
th2_54=ones(length(len))*max(x_54_);
hold on;plot(x_54,y_54_,'rs-');
title(['GBP Fluo 5.4, ',int2str(min(x_54_)),'-',int2str(max(x_54_))]);
find(deriv==0)

subplot(2,2,3);
hold on;plot(x_51,y_51_,'gd-');
deriv_51=diff(y_51_,order);
deriv1=diff(y_51_,1);
x_51_=x_51(1:end-2);
x_51_=x_51_(and(deriv1(1:end-1)>0.0099,abs(deriv_51)<=0.0099));
hold on;
th1_51=ones(length(len))*min(x_51_);
th2_51=ones(length(len))*max(x_51_);
stem(x_51(1:end-order),deriv_51,'gd-');
title(['GBP Fluo 5.1, ',int2str(min(x_51_)),'-',int2str(max(x_51_))]);
find(deriv==0)

saveas(f2,'Figures/GFPCurvesSeparate.bmp','bmp');

f1=figure(1);
hold on;
plot(x,y_,x_54,y_54_,'rs-',x_51,y_51_,'gd-')
plot(th1,len,th2,len,th1_54,len,'rs-',th2_54,len,'rs-',th1_51,len,'gd-',th2_51,len,'gd-');
legend(['GBP Fluo 5.6 linear:',int2str(min(x_56_)),'-',int2str(max(x_56_))],['GBP Fluo 5.4 linear:',int2str(min(x_54_)),'-',int2str(max(x_54_))],['GBP Fluo 5.1 linear:',int2str(min(x_51_)),'-',int2str(max(x_51_))] );
hold on;
stem(x(1:end-order),deriv);
stem(x_54(1:end-order),deriv_54,'rs-');
stem(x_51(1:end-order),deriv_51,'gd-');

saveas(f1,'Figures/GFPCurvesAllInOne.bmp','bmp');