function [Wsum]=SumWeigthsPerPosition(Pos,W,s,e)

Wsum=zeros(e-s,1);

count=1;
for i=s:1:e
    Wsum(count)=sum(abs(W(Pos==i)));
    count=count+1;
end

end