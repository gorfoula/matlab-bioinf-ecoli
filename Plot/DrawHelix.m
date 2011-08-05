function [text_handle]=DrawHelix(Seq,index,pos)

r=1;
aas=1:length(Seq);
degrees=(aas-1)*100;
circleDegr=mod(degrees,360);
radius=circleDegr.*(pi/180);
xlim([-2 2]);
ylim([-2 2]);
step=.5;
[AAs SelSeq SelSeqend]=AARepresentation('int',1,length(Seq),{Seq}); % assign amino acid IDs
[PolScale SumPol]=Polarity(AAs,'Z');
text_handle=zeros(1,length(Seq));
for i=1:length(Seq)
    x=r.*cos(radius);
    y=r.*sin(radius);
    if(exist('index','var') & sum(i==index)>0)
        switch Seq(i)
            case {'K','H','R'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','b','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','EdgeColor','red');
            case {'D','E'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','EdgeColor','red');
            case {'I','V','L','F','C','M','A'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','y','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','EdgeColor','red');
            otherwise
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','k','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','EdgeColor','red');
        end
    else
        switch Seq(i)
            case {'K','H','R'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','b','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
            case {'D','E'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','r','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
            case {'I','V','L','F','C','M','A'}
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','y','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
            otherwise
                th=text(x(i),y(i),[Seq(i),'(',num2str(i),')',num2str(PolScale(i))],'Color','k','HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
        end
    end
    if(mod(i,18)==0)
        r=r+step;
        xlim([-r-step r+step]);
        ylim([-r-step r+step]);
    end
text_handle(i)=th;
% pause;
end
% title([num2str(Profile_Z_i),' | ',num2str(Profile_Z_o)])

end