function [Not_FOUND_m2z Not_FOUND_byfeatures Count_OnePep_OneCharge]=Plot3D(EXP,THEOR,NONDET,maxPep)
close all;
ALL=cell(3,1);
% ID | Sequence |LEN|MASS|MAX CHARGE|m2z|GRAVY|Polarity|SerialNum
% MASS|GRAVY|MAX CHARGE|SerialNum
ALL{1,1}=CellTable2Double(EXP(:,[4 7 5 6 9]));
ALL{2,1}=CellTable2Double(THEOR(:,[4 7 5 6 9]));
ALL{3,1}=CellTable2Double(NONDET(:,[4 7 5 6 9]));
texts_={'Experimental' 'Theoretical' 'Non Detected'};
%%%%%
temp=hot;
hot_inv=temp(end:-1:1,:);
y_=(0:.1:1)*max(ALL{1,1}(:,1));y_ctrs=mean([y_(2:end);y_(1:end-1)]);
x_=-1:(2/(length(y_)-1)):1;x_ctrs=mean([x_(2:end);x_(1:end-1)]);
FREQ=zeros(size(x_,2),size(x_,2),3);
NumOfDetPept=1:maxPep;
PepDetectableArea=zeros(2,maxPep);
bins=20;
min_m2z=min([ALL{1,1}(:,4);ALL{2,1}(:,4);ALL{3,1}(:,4)]);
max_m2z=max([ALL{1,1}(:,4);ALL{2,1}(:,4);ALL{3,1}(:,4)]);
m2z_x=min_m2z:(max_m2z-min_m2z)./(bins-1):max_m2z;
mass2charge=zeros(2,length(m2z_x));
for i=1:3
    dat = ALL{i};
    PROTEINS=ALL{i,1}(end,end);
    n=zeros(length(x_),length(y_));
    Count_Non_m2z=0;    Count_OnePep_OneCharge=0;
    Not_FOUND_m2z=zeros(PROTEINS,1); NumDetPeptides=zeros(1,PROTEINS); NumNonDetPeptides=zeros(1,PROTEINS);
    for j=1:PROTEINS %% for each protein
        sub_data=dat(dat(:,end)==j,1:3); %% choose peptides of protein j
        [n_ x]= hist3(sub_data(:,1:2),'Ctrs',{y_;x_}); % 2D hist data;
        if(i==2 || i==3)  %% only for theoretical and non detected list of peptides
            detectablePeptides=(n_>0 & FREQ(:,:,1)>0);
            NotdetectablePeptides=(n_>0 & FREQ(:,:,1)==0);
            NumNonDetPeptides(j)=sum(sum(NotdetectablePeptides));
            NumDetPeptides(j)=sum(sum(detectablePeptides));
        end    
        if(sum(sum(n_))>0)
            n=n_+n; %% Add probability of a peptide being in a specific area
%             n=(n_./sum(sum(n_)))+n; %% Add probability of a peptide being in a specific area
            if(sum(sum(n_))==1)
                Count_OnePep_OneCharge=Count_OnePep_OneCharge+sum(sub_data(:,3)==1);
            end
        else
            Count_Non_m2z=Count_Non_m2z+1;
            Not_FOUND_m2z(j)=1;
        end    
    end    
    n=(n./size(dat,1))*100;
    FREQ(:,:,i)=n;
    
    figure(1);subplot(1,3,i);
    box off;
    xlabel('GRAVY');
    ylabel('MASS');
    colormap(jet);

    h=pcolor(x{1,2},x{1,1},n);
    set(h, 'zdata', ones(size(n)) *(i-1)*25);
    text(1.1,15,((i-1)*25)+5,[texts_{i}],'Color','k','FontWeight','bold');
    colorbar;
    Not_FOUND_byfeatures=sum(NumNonDetPeptides>0);
    display(['Number of tryptic Peptides not Found: ',num2str(sum(NumNonDetPeptides))]);
    display(['Number of Protein / all Tryptic Peptides not detected: ',num2str(Not_FOUND_byfeatures)]);
    
    NonDetectableArea=FREQ(:,:,1)<=0;
    Y_square=y_'*ones(size(y_));
    X_square=ones(size(x_))'*x_;
    [row col]=find(NonDetectableArea(1:end-1,1:end-1)==1);
    for sq=1:1:length(row)
        X_cur=[X_square(row(sq),col(sq)) X_square(row(sq),col(sq)) X_square(row(sq),col(sq)+1) X_square(row(sq),col(sq)+1)];
        Y_cur=[Y_square(row(sq),col(sq)) Y_square(row(sq)+1,col(sq)) Y_square(row(sq)+1,col(sq)) Y_square(row(sq),col(sq))];
        Z_cur=[(i-1)*25 (i-1)*25 (i-1)*25 (i-1)*25];
        patch(X_cur,Y_cur,Z_cur,Z_cur,'EdgeColor','r','LineWidth',1,'FaceColor','w');
    end

    if(i==2 || i==3)  %% only for theoretical and non detected list of peptides     
        mass2charge(i-1,:)=hist(dat(:,4),m2z_x)./length(dat(:,4));
        PepDetectableArea(i-1,:)=hist(NumDetPeptides,NumOfDetPept)./length(NumDetPeptides);
        PepDetectableArea_mn(i-1)=mean(NumDetPeptides);
        PepDetectableArea_mx(i-1)=max(PepDetectableArea(i-1,:));
        
        [func{i},xaxis{i}] = ecdf(NumDetPeptides);
    end
end
%%%%
figure;
[h6]=FigureLegends(m2z_x,mass2charge',6,'Mass 2 Charge','Peptides',[],{'Detected' 'NonDetected'},'b',{'-','';':','o'});

figure(4);plot(xaxis{2},func{2},'r',xaxis{3},func{3},'c')
legend('Detected','NonDetected');
xlabel('Peptides');
title('Cumulative Distribution of Number of Detectable Peptides');
figure;
[h5]=FigureLegends(NumOfDetPept,PepDetectableArea',5,'Number of Detectable Peptides','Proteins',[],{'Detected' 'NonDetected'},'b',{'-','';':','o'});

hold on;
y_axis=(0:.1:1)*PepDetectableArea_mx(1);
x_axis=ones(length(y_axis),1);
plot(x_axis*PepDetectableArea_mn(1),y_axis,'b.-');
text(PepDetectableArea_mn(1),PepDetectableArea_mx(1),num2str(PepDetectableArea_mn(1)),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');
y_axis=(0:.1:1)*PepDetectableArea_mx(2);
x_axis=ones(length(y_axis),1);
plot(x_axis*PepDetectableArea_mn(2),y_axis,'r.-');
text(PepDetectableArea_mn(2),PepDetectableArea_mx(2),num2str(PepDetectableArea_mn(2)),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');
hold off;



figure;SimpleHist(DETECTED,NOTDETECTED,[{'Detected'},{'NonDetected'}]);
ylabel('Number of Detectable Peptides');
xlabel('Peptides');
title('Cumulative Distribution of Number of Detectable Peptides');


% [freq x]=hist3(MOLWexp(:,2:3),[bins,bins]);
% X=MultipleCopiesVector(x{1,2},bins);
% Y=MultipleCopiesVector(x{1,1},bins);
% figure(2);plot3(X,Y',freq);
% figure(3);h=surfc(x{1,2},x{1,1},freq);%% Mass,kD,freq
% axis([-1 0.8 0 50 -800 800]);
% figure(4);hist3(MOLWexp(:,2:3),[bins,bins]);
% 
% [freq x]=hist3(MOLWtheo(:,2:3),[bins,bins]);
% X=MultipleCopiesVector(x{1,2},bins);
% Y=MultipleCopiesVector(x{1,1},bins);
% figure(2);hold on;plot3(X,Y',freq*-1);hold off;
% figure(3);hold on;h=surfc(x{1,2},x{1,1},freq*-1);hold off; %% Mass,kD,freq
% axis([-1 0.8 0 50 -800 800]);
% figure(4);hold on;hist3(MOLWtheo(:,2:3),[bins,bins],'r');hold off;



end

function [X2D]=MultipleCopiesVector(X,copies)
X2D=zeros(copies,length(X));
while (copies>0)
    X2D(copies,:)=X;
    copies=copies-1;
end

end

function []=SimpleHist(Fn,Nfn,leg)

Fn(Fn==0)=10^(-23);
Nfn(Nfn==0)=10^(-23);

if(length(unique(Fn))==1)
    [f_f x_f]=hist(Fn);f_f=f_f./sum(f_f);
else
    x_f=min(Fn):abs(max(Fn)-min(Fn))/19:max(Fn);
    [f_f x_f]=hist(Fn);f_f=f_f./sum(f_f);
end
if(length(unique(Nfn))==1)
    [f_n x_nf]=hist(Nfn);f_n=f_n./sum(f_n);
else
    x_nf=min(Nfn):abs(max(Nfn)-min(Nfn))/19:max(Nfn);
    [f_n x_nf]=hist(Nfn);f_n=f_n./sum(f_n);
end
hold on;

bar(x_f,f_f,0.9,'b');
bar(x_nf,f_n,0.7,'r');

mn_Fn=mean(Fn);mx_Fn=max(f_f);
mn_Nfn=mean(Nfn);mx_Nfn=max(f_n);

y_axis=(0:.1:1)*mx_Fn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Fn,y_axis,'b.-');
text(mn_Fn,mx_Fn,num2str(mn_Fn),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');
y_axis=(0:.1:1)*mx_Nfn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Nfn,y_axis,'r.-');
text(mn_Nfn,mx_Nfn,num2str(mn_Nfn),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');

hold off;
ylabel('Percent of CEP Proteins (%)');
legend(leg,'Location','Best');

[p,h,s] = kruskalwallis([Fn;Nfn],[ones(size(Fn,1),1);ones(size(Nfn,1),1)],'off');     display(p);
title_=['Kruskal Wallis p-value: ',num2str(p)];
title(title_);

end