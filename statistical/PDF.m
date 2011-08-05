function [freq x mx mn]=PDF(Y,groups,x)

RVs=max(groups)+1;
freq=zeros(length(x),RVs);
mn=zeros(RVs,1);
mx=zeros(RVs,1);

for rv=1:1:RVs
    [freq(:,rv)]=hist(Y(groups==(rv-1)),x);
    freq(:,rv)=freq(:,rv)*100./sum(freq(:,rv));
    mx(rv)=max(freq(:,rv));
    mn(rv)=mean(Y(groups==(rv-1)));
end

end