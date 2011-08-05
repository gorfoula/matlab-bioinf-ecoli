function []=LenDist(IFile,leg)
close all;
[FileNames] = textread(IFile,'%s',-1,'delimiter','\t');
FILES=length(FileNames);
SPLen_all=[];MatLen_all=[];

[colors cols]=size(colormap);
col_index=30:floor(colors/(FILES-1))-1:colors;

for  i=1:1:FILES
    [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(FileNames{i},'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
    if(SPLen(1)>0)
        SPLen_all=[SPLen_all;SPLen];
    end
    MatLen_all=[MatLen_all;(TotLen-SPLen)];
    
end
x=min(SPLen_all):4:max(SPLen_all);
x2=min(MatLen_all):150:1100;

SPLen_all=[];MatLen_all=[];
for  i=1:1:FILES
    [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(FileNames{i},'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
    
    MatLen=TotLen-SPLen;
    
    if(SPLen(1)>0)
        [freq x]=hist(SPLen,x);
        freq=freq./length(SPLen);
        SPLen_all=[SPLen_all freq'];
    end
    
    [freq x2]=hist(MatLen,x2);
    freq=freq./length(MatLen);
    MatLen_all=[MatLen_all freq'];

end
BarHist(1,1,x,(SPLen_all*100),'Signal Peptide Length',leg,col_index);
% BarHist(2,2,x2,(MatLen_all*100),'Mature Domain Length',[leg 'Cytoplasmic'],col_index);
end

function []=BarHist(fig,dim,x,Len,titl,leg,col_index)

figure(fig);
colormap gray;
if (dim==3)
    h = bar3(x,Len);
    set(h,'CDataMapping','direct');
    AssignColorIndex(h,col_index);
elseif (dim==2)
    grey_colors=colormap;
    h = bar(x,Len);
else
    [row col]=size(Len);
    strength=1;
    for i=1:1:col
        grey_colors=colormap;
        h = bar(x,Len(:,i),strength,'FaceColor',grey_colors(col_index(i),:));hold on;
        strength=strength-0.2;
    end
    hold off;
end

if(exist('legend','var')==0)
    l=legend(leg,'Location','Best');
end
    title(titl);
    ylabel('Length');zlabel('Percent over all');
    if(length(leg)==3)
        ylim([10 55]);
    end

end

function []=AssignColorIndex(h,col_index)
Zdata = get(h(1),'Zdata');
cdata = get(h(1),'Cdata');
cdata_size=size(cdata);
bar_images=length(h);

cdata=zeros(cdata_size,bar_images);
for i=1:1:bar_images
    cdata=ones(cdata_size)*col_index(i);
    set(h(i),'Cdata',cdata);
end

end