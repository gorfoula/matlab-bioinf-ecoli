close all;
window=1000;

b=randn(10^6,1);
[WeightFactor x]=hist(b,window);
mask=WeightFactor./10^6;
plot(x,mask);

figure(1);
x=0:0.1:8;
y=gaussmf(x,[1 4]);
plot(x,y);