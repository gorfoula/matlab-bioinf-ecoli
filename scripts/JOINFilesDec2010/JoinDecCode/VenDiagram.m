function []=VenDiagram(ZeroOneTable,labels)
[samples Bits]=size(ZeroOneTable);

if(length(labels)~=Bits)
    display('Number of lables id not the same as classes');
    return;
end
Combintions=(2^Bits)-1;
BinRepresentation=zeros(Combintions,Bits);
AREAS=zeros(Combintions,1);
for c=1:Combintions
    BinRepresentation(c,:)=dec2binvec(c,Bits);
    Propagate=ones(1,samples)'*dec2binvec(c,Bits);
    AREAS(c)=sum(sum(not(xor(ZeroOneTable,Propagate)),2)==Bits);
end
A=zeros(1,Bits);
for a=1:Bits
    A(a)=sum( AREAS( BinRepresentation(:,a)==1 ) );
end

InterSectAreas=sum(BinRepresentation,2)>1;
I=AREAS(InterSectAreas)';

figure, axis equal, axis off
%Using the same areas as above, display the error optimization at each iteration. Get the output structure.
F = struct('Display', 'iter');
[H,S] = venn(A,I,F,'ErrMinMode','TotalError','FaceAlpha', 0.6);
%Now label each zone 
Values=[AREAS(not(InterSectAreas))' I];
for i = Bits+1:Combintions
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(Values(i)),'Color','k','FontWeight','bold');
end

for b = 1:Bits
    text(S.Position(b,1), S.Position(b,2)+S.Radius(b), [labels{b},10,num2str(A(b))],'Color','k','FontWeight','bold','VerticalAlignment','bottom');
    yes=IsInOtherTable(S,Bits,b);
    while(yes==1)
        S.ZoneCentroid(b,1)=S.ZoneCentroid(b,1)-5;
%         S.ZoneCentroid(b,2)=S.ZoneCentroid(b,2)+.5;
        yes=IsInOtherTable(S,Bits,b);
    end
    text(S.ZoneCentroid(b,1), S.ZoneCentroid(b,2), num2str(Values(b)),'Color','k','FontWeight','bold');
end

end

function [yes]=IsInOtherTable(S,Bits,b)
    yes=0;
    for cur_circ=1:1:Bits
        if(cur_circ~=b)
            yes=yes | IsInCircle(S.ZoneCentroid(b,:),S.Position(cur_circ,:),S.Radius(cur_circ));
        end
    end

end
function [yes]=IsInCircle(XY,c,r)

yes=(c(1)-XY(1))^2+(c(2)-XY(2))^2<=r^2;

end