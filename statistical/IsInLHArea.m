function [LH Y Y_]=IsInLHArea(x,S,a,b)

X=S(:,1);
Y=S(:,2);
Y_=a*exp(-(X-x(1))/(x(end)-x(1)))+b;
LH=Y>Y_ & X>=1000;
end