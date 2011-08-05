function []=PropertiesOfTrypticPeptides()
clear all;
close all;
MAX=20;
step=100;
col=[5 6];
% col=[5 6 8];
r=.008; %for my scatter function size of spots
x=400:step:1800;
% x=800:step:5200;
x_gravy=-3:.1:3;
% [Experimental_Peptides]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_Matched.txt','\n');
% [Properties]=CellTable2Double(Experimental_Peptides(2:end,6:7));
% [m2z PeptMass PossibleCharges]=Mass2Charge(Experimental_Peptides(2:end,5));
% m2z=Properties(:,1)./PossibleCharges;
% FileWriteTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPepti
% des_Matched_additional.txt',[Experimental_Peptides [[{'m/z'} {'Charges'}];Double2CellTable([m2z PossibleCharges])]],[],'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/CEP_NonDetected_Common(1)F.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));A=A(2:end,:);

[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));B=B(2:end,:);

[C]=ReadTable('JoinFiles/CompareTrypticPeptides/BL21_CEP_Common(6)F.txt_instr.txt','\n');
[C_Double]=CellTable2Double(C(2:end,col));C=C(2:end,:);

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Undetected'} {'detected'}];

Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;hold on;

%%%%% Peptides in LH area  %%%%%%%
y0=1.7;y1=.85;a=(y1-y0)/(exp(-1)-1);b=y0-a;
[LH_A]=IsInLHArea(x,A_Double(:,1:2),a,b);
[LH_B]=IsInLHArea(x,B_Double(:,1:2),a,b);
[LH_C]=IsInLHArea(x,C_Double(:,1:2),a,b);

Myscatter([A_Double(:,1);B_Double(:,1)],[A_Double(:,2);B_Double(:,2)],groups,r);
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%% DRAW LH AREA  %%%%%%%%%%%%%
y=a*exp(-(x-x(1))/(x(end)-x(1)))+b;
y_1000=a*exp(-(1000-x(1))/(x(end)-x(1)))+b;
y_1800=a*exp(-(1800-x(1))/(x(end)-x(1)))+b;
p=patch([x(x>=1000) 1800 1000 1000],[y(x>=1000) 3 3 y_1000],[zeros(1,length(x(x>=1000))) 0 0 0],'EdgeColor','g','FaceColor','g','FaceAlpha',0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off;
display('Undetected: ');length(unique(A(LH_A,1)))
display('detected: ');length(unique(B(LH_B,1)))
%%%%%%%%%%%%%  LH COVERAGE  %%%%%%%%%%%%%%%%%%%
[LH_coverage_A Index_A]=LHCoverage(A(:,1),LH_A);
[LH_coverage_B Index_B]=LHCoverage(B(:,1),LH_B);
%--plot
Y=[LH_coverage_A;LH_coverage_B];
groupsP=[zeros(length(LH_coverage_A),1);ones(length(LH_coverage_B),1)]; % this group is different from the above in size (we have proteins, not peptides)
LabeL=[{''} {'LH coverage (%)'} {'Proteins (%)'}]; %% title|xlabel|ylabel
FigureLegends(0:10:100,Y,groupsP,LabeL,LegenD,'b',{'-','o';':','d'},0);

proteins=length(LH_coverage_A);
display(['Undetected|Proteins tryptic peptides above line: ',num2str(sum(LH_coverage_A>0)),'|',num2str(sum(LH_coverage_A>0)*100./proteins),'%']);
display(['Undetected|Proteins with all tryptic peptides above line: ',num2str(sum(LH_coverage_A==100)),'|',num2str(sum(LH_coverage_A==100)*100./proteins),'%']);

proteins=length(LH_coverage_B);
display(['Undetected|Proteins tryptic peptides above line: ',num2str(sum(LH_coverage_B>0)),'|',num2str(sum(LH_coverage_B>0)*100./proteins),'%']);
display(['Undetected|Proteins with all tryptic peptides above line: ',num2str(sum(LH_coverage_B==100)),'|',num2str(sum(LH_coverage_B==100)*100./proteins),'%']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HEAT MAP  %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TABLE=[A_Double(:,1:2);B_Double(:,1:2)];
Names=[A(:,1);B(:,1)];
ids=[Index_A;Index_B];
MyHeatMap(TABLE,groups,ids,LegenD,Names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Write Table S11   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxID=max(C_Double(:,end));
INDEX_LH=zeros(maxID,1);INDEX_nLH=zeros(maxID,1);
id_unique_LH=unique(C_Double(LH_C,end));
id_unique_nLH=unique(C_Double(not(LH_C),end));
INDEX_LH(id_unique_LH)=1;
INDEX_nLH(id_unique_nLH)=1;
INDEX_both=find(INDEX_LH==1 & INDEX_nLH==1);
INDEX_haonly=find(INDEX_LH==1 & INDEX_nLH==0);
FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[C(1,:) {'Is large and Hydrophobic'} {'Detected in this study'}],[],'w');
for i=INDEX_both'
    cur_haindex=(C_Double(:,end)==i) & (LH_C==1);
    cur_notindex=(C_Double(:,end)==i) & (LH_C==0);
    label=cell(sum(cur_haindex),2);label(:,1)={'yes'};
    cur_table=C(logical([0;cur_haindex]),:);
    if(sum(strcmp(A(1:end,1),cur_table(1,1)))==0)
        label(:,2)={'yes'};
    end
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[cur_table label],[],'a');
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[C(logical([0;cur_notindex]),:) cell(sum(cur_notindex),1)],[],'a');
end
FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[C(1,:) {'Is large and Hydrophobic'}],[],'w');
for i=INDEX_haonly'
    cur_haindex=(C_Double(:,end)==i) & (LH_C==1);
    cur_notindex=(C_Double(:,end)==i) & (LH_C==0);
    label=cell(sum(cur_haindex),2);label(:,1)={'yes'};
    cur_table=C(logical([0;cur_haindex]),:);
    if(sum(strcmp(B(1:end,1),cur_table(1,1)))>0)
        label(:,2)={'yes'};
    end
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[cur_table label],[],'a');
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[C(logical([0;cur_notindex]),:) cell(sum(cur_notindex),1)],[],'a');
end

[P1 P2 P3 P4 P5 FILES]=FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[1 1],[1 1],[6]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt',FILES{1},[1 2],[1 1],[]);
FishLines('JoinFiles/k12/uniprot-organism_83333_SubLocFixed.txt','JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[1 1],[1 1]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/CompareTrypticPeptides/(DB)uniprot-organism_83333_SubLocFixed(for)CEP_Large_Pho_Peptides_only(16)F.txt',[1 6],[1 1]);

clear P1 P2 P3 P4 P5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A]=ReadTable('JoinFiles/CompareTrypticPeptides/Non_DetByUs_CEP(1)F.txt_count.txt','\n');
% [A_Double]=CellTable2Double(A(2:end,col));
% [B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
% [B_Double]=CellTable2Double(B(2:end,col));

[A]=ReadTable('JoinFiles/CompareTrypticPeptides/Non_DetByUs_CEP(1)F.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

[C]=ReadTable('JoinFiles/CompareTrypticPeptides/BL21_CEP.txt_count.txt','\n');
[C_Double]=CellTable2Double(C(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Undetected'} {'Detected'}];
% Colors=[0 0 0;0 1 0];
Colors=[0 0 0;1 0 1];

norm_factor=max(C_Double);

C_Double_norm=[C_Double(:,1)./norm_factor(1) C_Double(:,2)./norm_factor(2)];
B_Double_norm=[B_Double(:,1)./norm_factor(1) B_Double(:,2)./norm_factor(2)];
A_Double_norm=[A_Double(:,1)./norm_factor(1) A_Double(:,2)./norm_factor(2)];


[idx3,ctrs,sumd,D] = kmeans(C_Double_norm,2,'distance','sqEuclidean','Replicates',100);ctrs_=[ctrs(:,1)*norm_factor(1) ctrs(:,2)*norm_factor(2)];
figure;hold on;
scatter(C_Double(idx3==1,1),C_Double(idx3==1,2),'MarkerEdgeColor','k','Marker','o','SizeData',72^(1/2));
scatter(C_Double(idx3==2,1),C_Double(idx3==2,2),'MarkerEdgeColor','k','Marker','.','SizeData',72);

[indx_A]=Classify_Kmeans(A_Double_norm,ctrs);ha_rect_A=(A_Double(:,1)>1000 & A_Double(:,2)>1);
[indx_B]=Classify_Kmeans(B_Double_norm,ctrs);ha_rect_B=(B_Double(:,1)>1000 & B_Double(:,2)>1);
[m ha]=max(ctrs(:,1));ha_rect=(C_Double(:,1)>1000 & C_Double(:,2)>1);
scatter(A_Double(indx_A==ha,1),A_Double(indx_A==ha,2),'MarkerEdgeColor','g','Marker','x','SizeData',72);
% scatter(A_Double(indx_A==ha_rect_A,1),A_Double(indx_A==ha_rect_A,2),'MarkerEdgeColor',Colors(2,:),'Marker','x','SizeData',72);
xlabel('Mean m/z');ylabel('Mean GRAVY');
% plot(ctrs_(:,1),ctrs_(:,2),'r.','LineWidth',2);

display(['Non Detected: ',num2str(round(sum(indx_A==ha_rect_A)*100./length(indx_A))),'%']);
display(['Detected: ',num2str(round(sum(indx_B==ha_rect_B)*100./length(indx_B))),'%']);

ha_rect=(B_Double(:,1)>=1000 & B_Double(:,2)>=1);
maxID=max(B_Double(:,end));
INDEX_ha=zeros(maxID,1);
INDEX_nha=zeros(maxID,1);
id_unique_ha=unique(B_Double(ha_rect,end));
id_unique_nha=unique(B_Double(not(ha_rect),end));
INDEX_ha(id_unique_ha)=1;
INDEX_nha(id_unique_nha)=1;
INDEX_both=find(INDEX_ha==1 & INDEX_nha==1);
INDEX_haonly=find(INDEX_ha==1 & INDEX_nha==0);

length(INDEX_both)./length(unique(B_Double(:,end)))

FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[C(1,:) {'ha'}],[],'w');
for i=INDEX_both'
    cur_haindex=(C_Double(:,end)==i) & (ha_rect==1);
    cur_notindex=(C_Double(:,end)==i) & (ha_rect==0);
    label=cell(sum(cur_haindex),1);label(:,:)={'ha'};
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[C(logical([0;cur_haindex]),:) label],[],'a');
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[C(logical([0;cur_notindex]),:) cell(sum(cur_notindex),1)],[],'a');
end
FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[C(1,:) {'ha'}],[],'w');
for i=INDEX_haonly'
    cur_haindex=(C_Double(:,end)==i) & (ha_rect==1);
    cur_notindex=(C_Double(:,end)==i) & (ha_rect==0);
    label=cell(sum(cur_haindex),1);label(:,:)={'ha'};
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[C(logical([0;cur_haindex]),:) label],[],'a');
    FileWriteTable('JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[C(logical([0;cur_notindex]),:) cell(sum(cur_notindex),1)],[],'a');
end

FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides.txt',[2 1],[1 1]);
FishLines('JoinFiles/BL2119/TableS4_BL21_DE3_ProteomeV8.txt','JoinFiles/CompareTrypticPeptides/CEP_Large_Pho_Peptides_only.txt',[2 1],[1 1]);
% FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/CompareTrypticPeptides/Non_DetByUs_CEP(1)F.txt_instr.txt',[2 1],[1 1]);
% FishLines('JoinFiles/CompareTrypticPeptides/(DB)TableS3_K12_ProteomeV13(for)Non_DetByUs_CEP(1)F.txt_instr(2)F.txt','JoinFiles/CompareTrypticPeptides/temp.txt',[1 1],[1 1]);
hold off;

f=figure;
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0,[f 2 1 1]);
Y=[A_Double(:,2);B_Double(:,2)];
FigureLegends(x_gravy,Y,groups,[{''} {'GRAVY'} {'Proteins (%)'}],LegenD,'p',{'-','o';':','d'},0,[f 2 1 2]);

TABLE=[A_Double(:,1:2);B_Double(:,1:2)];
Names=[A(2:end,1);B(2:end,1)];
ids=[CellTable2Double(A(2:end,end));CellTable2Double(B(2:end,end))];
[f1,n,x_,y_]=MyHeatMap(TABLE,groups,ids,LegenD,Names,Colors);

figure;hold on;

scatter(B_Double(:,1),B_Double(:,2),'MarkerEdgeColor',Colors(2,:),'Marker','o','SizeData',72^(1/2));
scatter(A_Double(:,1),A_Double(:,2),'MarkerEdgeColor',Colors(1,:),'Marker','.','SizeData',72);hold on;
xlabel('Mean m/z');ylabel('Mean GRAVY');
legend([LegenD(1,2);LegenD(1,1)]);

Levels=[];
Properties=cell(1,2);
Properties{1,1}=A_Double;
Properties{2,1}=B_Double;
for i=1:1:size(n,3);
    [X,Y] = meshgrid(x_,y_);
    [C,h]=contour(Y,X,n(:,:,i),25);
    cur_Levels=get(h,'LevelList');
    Levels=[Levels cur_Levels(1)];
    set(h,'LevelList',Levels);
    set(h,'LineWidth',1);
    set(h,'LineColor',Colors(i,:));
    
    [freq]=hist(Properties{i,1}(:,1),y_);freq=(freq./sum(freq))+min(x_);
%     plot(y_,freq,'Color',Colors(i,:));
    y_axis=x_;
    x_axis=ones(length(y_axis),1)*mean(Properties{i,1}(:,1));
    plot(x_axis,y_axis,'Color',Colors(i,:),'LineStyle',':');
    y_axis=y_;
    x_axis=ones(length(y_axis),1)*mean(Properties{i,1}(:,2));
    plot(y_axis,x_axis,'Color',Colors(i,:),'LineStyle',':');
end

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/BL21_Det_CEP.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/BL21_NonDet_CEP.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

x=800:step:5200;
LabeL=[{''} {'Mass'} {'Tryptic Peptides (%)'}]; %% title|xlabel|ylabel
LegenD=[{'BL21 Det'} {'BL21 Undetected'}];

figure;
scatter(A_Double(:,1),A_Double(:,2),'MarkerEdgeColor','r','Marker','.','SizeData',72);hold on;
scatter(B_Double(:,1),B_Double(:,2),'MarkerEdgeColor','b','Marker','o','SizeData',72^(1/2));
xlabel('Mass');ylabel('GRAVY');
legend(LegenD);

y0=2;
y1=-0.25;
a=(y1-y0)/(exp(-1)-1);
b=y0-a;
y=a*exp(-(x-x(1))/(x(end)-x(1)))+b;

plot(x,y,'k');

YB=B_Double(:,2);
XB=B_Double(:,1);
YB_=a*exp(-(XB-x(1))/(x(end)-x(1)))+b;
YA=A_Double(:,2);
XA=A_Double(:,1);
YA_=a*exp(-(XA-x(1))/(x(end)-x(1)))+b;
scatter(A_Double(YA>YA_,1),A_Double(YA>YA_,2),'MarkerEdgeColor','r','Marker','.','SizeData',72);hold on;
hold off;
display('BL21 CEP Detected: ');length(unique(A_Double(YA>YA_,3)))
display('BL21 CEP Non Detected: ');length(unique(B_Double(YB>YB_,3)))

index=YA>YA_;
count=0;total_count=0;
for i=unique(A_Double(:,3))'
    cur_prot_peptides=index(A_Double(:,3)==i);
    if(sum(cur_prot_peptides)==length(cur_prot_peptides))
        count=count+1;
    end
    if(sum(cur_prot_peptides)>0)
        total_count=total_count+1;
    end

end
display(['Detected|Proteins tryptic peptides above line: ',num2str(total_count),'/',num2str(length(unique(A_Double(:,3))))]);
display(['Detected|Proteins with all tryptic peptides above line: ',num2str(count),'/',num2str(length(unique(A_Double(:,3))))]);

index=YB>YB_;
count=0;total_count=0;
name_index=zeros(size(index,1),1);
print_table=zeros(size(index,1),2);
for i=unique(B_Double(:,3))'
    cur_prot_peptides=index(B_Double(:,3)==i);
    if(sum(cur_prot_peptides)==length(cur_prot_peptides))
        count=count+1;
    end
    if(sum(cur_prot_peptides)>0)
        total_count=total_count+1;
        name_index(find(B_Double(:,3)==i,1,'first'))=1;
        print_table(find(B_Double(:,3)==i,1,'first'),:)=[size(cur_prot_peptides,1) sum(cur_prot_peptides)];
    end
end
display(['NonDetected|Proteins tryptic peptides above line: ',num2str(total_count),'/',num2str(length(unique(B_Double(:,3))))]);
display(['NonDetected|Proteins with all tryptic peptides above line: ',num2str(count),'/',num2str(length(unique(B_Double(:,3))))]);

TABLE=[A_Double(YA>YA_,1:2);B_Double(YB>YB_,1:2)];
groups=[zeros(sum(YA>YA_),1);ones(sum(YB>YB_),1)];
Names=[A(logical([0;YA>YA_]),1);B(logical([0;YB>YB_]),1)];
Y=[A_Double(YA>YA_,1);B_Double(YB>YB_,1)];

print_table=[[{'Total detectable Tryptic Peptides'} {'Total Hydrophobic Peptides'}];Double2CellTable(print_table)];
FileWriteTable('JoinFiles/CompareTrypticPeptides/Proteins_Large_Pho_Peptides_NonDet_CEP.txt',[B(logical([1;name_index]),[1 3]) print_table(logical([1;name_index]),:)],[],'w');

% groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);

ids=[CellTable2Double(A(logical([0;YA>YA_]),end));CellTable2Double(B(logical([0;YB>YB_]),end))];
MyHeatMap(TABLE,groups,ids,LegenD,Names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/2ExperimentalTrypticPeptides_CEP.txt_count.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));

[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Experimental'} {'Theoretical'}];

Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);

TABLE=[A_Double(:,1:2);B_Double(:,1:2)];
Names=[A(2:end,1);B(2:end,1)];
ids=[CellTable2Double(A(2:end,end));CellTable2Double(B(2:end,end))];
[f1,n,x_,y_]=MyHeatMap(TABLE,groups,ids,LegenD,Names);


figure;hold on;

scatter(A_Double(:,1),A_Double(:,2),'MarkerEdgeColor','k','Marker','.','SizeData',72);hold on;
scatter(B_Double(:,1),B_Double(:,2),'MarkerEdgeColor','m','Marker','o','SizeData',72^(1/2));

xlabel('Mean m/z');ylabel('Mean GRAVY');
legend(LegenD);

Levels=[];
Colors=[0 0 0;1 0 1];
for i=1:1:size(n,3);
    [X,Y] = meshgrid(x_,y_);
    [C,h]=contour(Y,X,n(:,:,i),30);
    cur_Levels=get(h,'LevelList');
    Levels=[Levels cur_Levels(1)];
    set(h,'LevelList',Levels);
    set(h,'LineWidth',1);
    set(h,'LineColor',Colors(i,:));
end
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/CYTO(1)F.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/CEP(1)F.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'CEP'} {'CYTO'}];
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;
scatter(A_Double(:,1),A_Double(:,2),'MarkerEdgeColor','r','Marker','.','SizeData',72);hold on;
scatter(B_Double(:,1),B_Double(:,2),'MarkerEdgeColor','b','Marker','o','SizeData',72^(1/2));
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);

y0=1.8;
y1=.8;
a=(y1-y0)/(exp(-1)-1);
b=y0-a;
y=a*exp(-(x-x(1))/(x(end)-x(1)))+b;
plot(x,y,'k');

Y=A_Double(:,2);
X=A_Double(:,1);
IDENT=A_Double(:,1);
Y_=a*exp(-(X-x(1))/(x(end)-x(1)))+b;
scatter(A_Double(Y>Y_,1),A_Double(Y>Y_,2),'MarkerEdgeColor','r','Marker','.','SizeData',72);hold on;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_out.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));

[B]=ReadTable('JoinFiles/CompareTrypticPeptides/Non_DetByUs_CEP(1)F.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected'} {'Undetected'}];

Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;

scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2),'MarkerEdgeColor','k','Marker','.','SizeData',72);hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','m','Marker','.','SizeData',72);

xlabel('m/z');ylabel('GRAVY');
legend(LegenD);

TABLE=[A_Double(:,:);B_Double(:,:)];
Names=[A(2:end,:);B(2:end,:)];
ids=[CellTable2Double(A(2:end,end));CellTable2Double(B(2:end,end))];
OUT_OF_DETECTABLE_AREA=MyHeatMap(TABLE,groups,ids,LegenD,Names);
FileWriteTable('JoinFiles/CompareTrypticPeptides/OUT_OF_DETECTABLE_AREA.txt',OUT_OF_DETECTABLE_AREA,[],'w');
FishLines('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','JoinFiles/CompareTrypticPeptides/OUT_OF_DETECTABLE_AREA.txt',[2 1],[1 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_Cyto(1)F.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

[C]=ReadTable('JoinFiles/CompareTrypticPeptides/NotDetByUs_NotDetmRNALevel(1)F.txt_instr.txt','\n');
[C_Double]=CellTable2Double(C(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CYTO'} {'Detected CEP'} {'UnDetected CEP'}];
Y=[A_Double(:,1);B_Double(:,1);C_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1);ones(length(C_Double(:,1)),1)*2];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d';'.','*'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2),'MarkerEdgeColor','y','Marker','s','SizeData',8);hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','b','Marker','.','SizeData',25);
scatter(C_Double(C_Double(:,1)>400,1),C_Double(C_Double(:,1)>400,2),'MarkerEdgeColor','c','Marker','.','SizeData',25);hold off;
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_zero_mis.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_mis.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP experimental (0 misscleavages)'} {'Detected CEP experimental (1-3 misscleavages)'}];
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

[C]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_mis.txt','\n');
[C_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP theoretical'} {'Detected CEP experimental (cleaved completely)'} {'Detected CEP experimental (1-3 misscleavages'}];
Y=[A_Double(:,1);B_Double(:,1);C_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1);ones(length(C_Double(:,1)),1)*2];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d';'.','*'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(C_Double(C_Double(:,1)>400,1),C_Double(C_Double(:,1)>400,2),'MarkerEdgeColor','y','Marker','d');
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

[C]=ReadTable('JoinFiles/CompareTrypticPeptides/ExperimentalTrypticPeptides_CEP.txt_zero_mis.txt','\n');
[C_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'m/z'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP theoretical'} {'Detected CEP experimental (cleaved completely)'} {'Detected CEP experimental (0 misscleavages)'}];
Y=[A_Double(:,1);B_Double(:,1);C_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1);ones(length(C_Double(:,1)),1)*2];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d';'.','*'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(C_Double(C_Double(:,1)>400,1),C_Double(C_Double(:,1)>400,2),'MarkerEdgeColor','y','Marker','d');
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/NotDetByUs_mRNALevel(1)F.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));
LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'NonDetected CEP (but have trancript)'}];
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/NotDetByUs_NotDetmRNALevel(1)F.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'NonDetected CEP (also no trancript)'}];
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
xlabel('m/z');ylabel('GRAVY');
legend(LegenD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_instr.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_Cyto(1)F.txt_instr.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'Detected CYTO'}];
Y=[A_Double(:,1);B_Double(:,1)];
groups=[zeros(length(A_Double(:,1)),1);ones(length(B_Double(:,1)),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
figure;
scatter(A_Double(A_Double(:,1)>400,1),A_Double(A_Double(:,1)>400,2));hold on;
scatter(B_Double(B_Double(:,1)>400,1),B_Double(B_Double(:,1)>400,2),'MarkerEdgeColor','r','Marker','.');
legend(LegenD);
xlabel('m/z');ylabel('GRAVY');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/Theoretical_CHECK_IDs_Theo(2)F_CEP.txt_count.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP New'} {'NonDetected CEP old'}];
Y=[A_Double;B_Double];
groups=[zeros(length(A_Double),1);ones(length(B_Double),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/NotDetByUs_mRNALevel(1)F.txt_count.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'NonDetected CEP (but have trancript)'}];
Y=[A_Double;B_Double];
groups=[zeros(length(A_Double),1);ones(length(B_Double),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/NotDetByUs_NotDetmRNALevel(1)F.txt_count.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'NonDetected CEP (also no trancript)'}];
Y=[A_Double;B_Double];
groups=[zeros(length(A_Double),1);ones(length(B_Double),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_CEP(1)L.txt_count.txt','\n');
[A_Double]=CellTable2Double(A(2:end,col));
[B]=ReadTable('JoinFiles/CompareTrypticPeptides/DetByUs_Cyto(1)F.txt_count.txt','\n');
[B_Double]=CellTable2Double(B(2:end,col));

LabeL=[{''} {'Number of Detectable Tryptic Peptides'} {'Proteins (%)'}]; %% title|xlabel|ylabel
LegenD=[{'Detected CEP'} {'Detected CYTO'}];
Y=[A_Double;B_Double];
groups=[zeros(length(A_Double),1);ones(length(B_Double),1)];
FigureLegends(x,Y,groups,LabeL,LegenD,'p',{'-','o';':','d'},0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [LH_coverage Index_A]=LHCoverage(A,LH_A)
% LH_A: index of which rows (proteins) in A are in LH area

    [countTP_A Index_A]=Douplicates('..txt',[0 1],A);
    [countTP_LHA_cell]=Douplicates('..txt',[0 1],StrMatrix2CellMatrix(num2str(Index_A(LH_A))));
    Proteins_with_LH_index=unique(Index_A(LH_A));
    countTP_LHA=zeros(size(countTP_A,1),1);
    countTP_LHA(Proteins_with_LH_index)=CellTable2Double(countTP_LHA_cell(:,2));
    LH_coverage=(countTP_LHA./CellTable2Double(countTP_A(:,2)))*100;

end