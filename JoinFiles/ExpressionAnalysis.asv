function []=ExpressionAnalysis(if_ProteinLevel,if_mRNALevel,sel)
close all;

[TABLE_mRNA Quant_mRNA Check_mRNA category_mRNA EmptyCells_mRNA NotDetected_mRNA Detected_Us_mRNA USDETECTED_mRNA USNOTDETECTED_mRNA index_mRNA]=PlotFounVsNotFound(if_mRNALevel,sel,7);
[TABLE_P Quant_P Check_P category_P EmptyCells_P NotDetected_P Detected_Us_P USDETECTED_P USNOTDETECTED_P index]=PlotFounVsNotFound(if_ProteinLevel,sel,13);

OutofRange=Check_P(:,end-2);
Empty=EmptyCells_mRNA & not(NotDetected_P);
LowDet_mRNA=(sum(Quant_mRNA(:,1:end-1),2)<=1) & (Check_mRNA(:,6)==1) & (Quant_mRNA(:,end)<1);
NotDet_mRNA=(sum(Quant_mRNA(:,1:end-1),2)<=1) & not(EmptyCells_mRNA) & not(LowDet_mRNA);
Transcribed=(not(NotDet_mRNA) & not(EmptyCells_mRNA));
NotDet_P=NotDetected_P;
%%%%%%%%%%%%%  CEP  %%%%%%%%%%%%%%%%%%%%%%
CEP=(category_mRNA==1);
CEP_Empty=(EmptyCells_mRNA & category_mRNA==1);
CEP_notTrancribed=(NotDet_mRNA & category_mRNA==1);
CEP_LowTrancribed=(LowDet_mRNA & category_mRNA==1);
CEP_Trancribed=(Transcribed & category_mRNA==1);
CEP_Detected=(not(NotDet_P) & category_mRNA==1);
%%%%%%%%%%%%%%%%  US  CYTO %%%%%%%%%%%%%%%%%%%%%
CYTO=(category_mRNA==0);
Cyto_Us_Detected=(Detected_Us_mRNA & category_mRNA==0);
%%%%%%%%%%%%%%%%  US  Non Det %%%%%%%%%%%%%%%%%%%%%
CEP_Us_Not_Detected_Empty=(not(Detected_Us_mRNA) & category_mRNA==1 & CEP_Empty);
CEP_Us_Not_Detected=(not(Detected_Us_mRNA) & category_mRNA==1);
CEP_Us_Not_Detected_not_mRNA=CEP_Us_Not_Detected & CEP_notTrancribed;
CEP_Us_Not_Detected_low_mRNA=CEP_Us_Not_Detected & CEP_LowTrancribed;
CEP_Us_Not_Detected_empty=CEP_Us_Not_Detected & CEP_Empty;
CEP_Us_Not_Detected_trancribed=CEP_Us_Not_Detected & (CEP_Trancribed) & not(CEP_Empty) & not(CEP_LowTrancribed) & not(OutofRange);
%%%%%%%%%%%%%%%%  US  Det %%%%%%%%%%%%%%%%%%%%%
CEP_Us_Detected=(Detected_Us_mRNA & category_mRNA==1);
CEP_Us_Detected_not_mRNA=CEP_Us_Detected & CEP_notTrancribed;
CEP_Us_Detected_low_mRNA=CEP_Us_Detected & CEP_LowTrancribed;
CEP_Us_Detected_empty=CEP_Us_Detected & CEP_Empty;
CEP_Us_Detected_trancribed=CEP_Us_Detected & (CEP_Trancribed | CEP_LowTrancribed) & not(CEP_Empty);
%%%%%%%%%%%%%%%%%%%%%%% mRNA abundance
prot_molecules_cell=Quant_mRNA(:,[7]);
index_molecules_cell=index_mRNA(:,[7]);
mRNA_molecules_cell=sum(prot_molecules_cell.*index_molecules_cell,2)./sum(index_molecules_cell,2);
isNAN=isnan(mRNA_molecules_cell);

det_=sum(Check_mRNA(:,[1 2 3]),2);
CEP_mRNA_Measured=(CEP_Trancribed & not(isNAN) & det_>0);
DETECTED=mRNA_molecules_cell(CEP_Us_Detected_trancribed & not(isNAN) & det_>0);
NOTDETECTED=mRNA_molecules_cell(CEP_Us_Not_Detected_trancribed & not(isNAN) & det_>0);

leg=[{['Detected by our study ','(',num2str(length(DETECTED)*100/sum(CEP_Us_Detected_trancribed)),'%)']} {['NOT Detected by our study','(',num2str(length(NOTDETECTED)*100/sum(CEP_Us_Not_Detected_trancribed)),'%)']}];
LabeL=[{'Cumulative Distribution of Number of mRNA molecules/cell'} {'mRNA molecules/cell'} {'Percent of CEP Proteins (%)'}];
SimpleHist(DETECTED,NOTDETECTED,LabeL,leg,[{[0:15:140]} {[0 15]}]);
%%%%%%%%%%%%%%%%%%%%%%% Protein mean abundance
prot_molecules_cell=Quant_P(:,[3 6 8]);
index_molecules_cell=index(:,[3 6 8]);
MEAN_molecules_cell=sum(prot_molecules_cell.*index_molecules_cell,2)./sum(index_molecules_cell,2);
isNAN=isnan(MEAN_molecules_cell);

det_=sum(Check_P(:,[1 2 3]),2);
CEP_Protein_Measured=(CEP_Detected & not(isNAN) & det_>0);
DETECTED=MEAN_molecules_cell(CEP_Us_Detected_trancribed & not(isNAN) & det_>0);
NOTDETECTED=MEAN_molecules_cell(CEP_Us_Not_Detected_trancribed & not(isNAN) & det_>0);

leg=[{['Detected by our study ','(',num2str(length(DETECTED)*100/sum(CEP_Us_Detected_trancribed)),'%)']} {['NOT Detected by our study','(',num2str(length(NOTDETECTED)*100/sum(CEP_Us_Not_Detected_trancribed)),'%)']}];
LabeL=[{'Cumulative Distribution of Number of protein molecules/cell'} {'protein molecules/cell'} {'Percent of CEP Proteins (%)'}];
SimpleHist(DETECTED,NOTDETECTED,LabeL,leg,[{[0:1000:12000]} {[0 6000]}]);
%%%%%%%%%%%%%%%%%%%  Correlation  %%%%%%%%%%%
CEP_Both=CEP_mRNA_Measured & CEP_Protein_Measured;
X=[mRNA_molecules_cell(CEP_Both) MEAN_molecules_cell(CEP_Both)];
figure;
scatter(X(:,1),X(:,2));
xlabel('mRNA molecules per cell');
ylabel('Protein molecules per cell');
RHO = corr(X);
title(['Pearson Correlation: ',num2str(RHO(1,2))]);
%%%%%%%%%%%%%%%%%%%%%%

display(['Cyto:  ',num2str(sum(CYTO))]);
display(['CEP:  ',num2str(sum(CEP))]);
display(['Empty:  ',num2str(sum(CEP_Empty))]);
display(['Not Transcribed: ',num2str(sum(CEP_notTrancribed))]);
display(['Low Transcribed: ',num2str(sum(CEP_LowTrancribed))]);
display(['Transcribed: ',num2str(sum(CEP_Trancribed))]);

Abudance=cell(size(TABLE_mRNA,1)-1,1);
Abudance(NotDet_mRNA,:)={'A'};
Abudance(LowDet_mRNA,:)={'L'};
Abudance(EmptyCells_mRNA,:)={'E'};
Abudance(not(EmptyCells_mRNA | LowDet_mRNA |NotDet_mRNA),:)={'P'};

FileWriteTable('JoinFiles/Transcriptomics/CEP.txt',TABLE_mRNA(logical([1;CEP]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/CYTO.txt',TABLE_mRNA(logical([1;CYTO]),:),[],'w');
% CEP_Us_Not_Detected_trancribed
FileWriteTable('JoinFiles/Transcriptomics/CEP_NonDetected.txt',TABLE_mRNA(logical([1;CEP_Us_Not_Detected]),:),[],'w');
FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/CEP_NonDetected.txt',[1 1],[1 1]);

FishLines('JoinFiles/k12/uniprot-organism_83333.txt','JoinFiles/Transcriptomics/NotDetByUs_Morethan10ProteinCopies.txt',[1 1],[1 0]);

FileWriteTable('JoinFiles/Transcriptomics/CEP_Trancribed.txt',TABLE_mRNA(logical([1;CEP_Trancribed]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/Non_DetByUs_CEP.txt',TABLE_mRNA(logical([1;CEP_Us_Not_Detected]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/DetByUs_Cyto.txt',TABLE_mRNA(logical([1;Cyto_Us_Detected]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/DetByUs_CEP.txt',TABLE_mRNA(logical([1;CEP_Us_Detected]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/NotDetByUs_mRNALevel.txt',[TABLE_mRNA(logical([1;CEP_Us_Not_Detected_trancribed]),:) TABLE_P(logical([1;CEP_Us_Not_Detected_trancribed]),:)],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/NotDetByUs_NotDetmRNALevel.txt',TABLE_mRNA(logical([1;CEP_Us_Not_Detected_not_mRNA]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/AllmRNAinfo.txt',[TABLE_mRNA [{'Abundance Level'};Abudance]],[],'w');
%[pvalue_X2 X2] = independenceTest(Check_Table, DESCR_MolPerCell' , [max(Check_Table) 2], 0);
end

function []=SimpleHist(Fn,Nfn,LabeL,leg,x)

bins=5;
tot=([Fn;Nfn]);

if(exist('x','var')==0)      %if expr is not defined
    x_tot=min(tot):abs(max(tot)-min(tot))/(bins -1):max(tot);
else
    x_tot=x{1};
    x_lim=x{2};
end

if(length(unique(Fn))==1)
    [f_f x_tot]=hist(Fn);f_f=f_f*100./sum(f_f);
else
    [f_f x_tot]=hist(Fn,x_tot);f_f=f_f*100./sum(f_f);
end
if(length(unique(Nfn))==1)
    [f_n x_tot]=hist(Nfn);f_n=f_n*100./sum(f_n);
else
    [f_n x_tot]=hist(Nfn,x_tot);f_n=f_n*100./sum(f_n);
end
groups=[zeros(length(Fn),1);ones(length(Nfn),1)];
Colors=FigureLegends(x_tot,[Fn;Nfn],groups,[{''} LabeL(2:end)],leg,'p',{'-','o';'--','d'},0,[2,1,1]);

hold on;
subplot(2,1,2);
hold on;
[func{1},xaxis{1}] = ecdf(Fn);
[func{2},xaxis{2}] = ecdf(Nfn);
semilogx(xaxis{1},func{1}*100,'Color',Colors(1,:),'MarkerFaceColor',Colors(1,:),'Marker','o','MarkerSize',2);
semilogx(xaxis{2},func{2}*100,'Color',Colors(2,:),'MarkerFaceColor',Colors(2,:),'Marker','d','MarkerSize',2);
% xlim([0 6000]);
xlim(x_lim);
title(LabeL{1});
xlabel(LabeL{2});
ylabel(LabeL{3});
legend(leg);
hold off;

end

function [TABLE Quant Check_Table category EmptyCells NotDetected Detected_Us USDETECTED USNOTDETECTED index]=PlotFounVsNotFound(FilePath,sel,plot_)

%%% sel:    0 for comparing found not found
%%%         1 rerun (do not read previous saved CEP files) compare found VS not found

if(sel==1)
    [TABLE category files COLS]=ExludeCyto(FilePath,1,1);
else
    [TABLE category files COLS]=ExludeCyto(FilePath,1,0);
end

[Quant index]=CellTable2Double(TABLE(2:end,files+COLS(1,1)+3:end));
[Check_Table]=CellTable2Double(TABLE(2:end,2:files+2));
Header=TABLE(1,files+COLS(1,1)+3:end);
Header_files=TABLE(1,2:files+1);
Detected_Us=Check_Table(:,files)==1;
EmptyCells=(sum(index,2)==0);
NotDetected=Check_Table(:,end)==0;
%  (DETECTED & (Measured At Zero | Not Measured)) | NOT DETECTED
Detected_Us_square = Detected_Us*ones(size(Detected_Us'));
USDETECTED=logical(Detected_Us_square(:,1:size(index,2)).*index);
USNOTDETECTED=logical(not(Detected_Us_square(:,1:size(index,2))).*index);

% for i=1:plot_
%     f=figure;
%     leg=[{[Header_files{end},'(',num2str(sum(USDETECTED(:,i))),')']};{['[NOT] ',Header_files{end},'(',num2str(sum(USNOTDETECTED(:,i))),')']}];
%     subplot(2,1,1);PlotQuant(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),leg);xlabel([Header{i},' (semilog axis)']);
%     subplot(2,1,2);SimpleHist(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),leg);xlabel(Header{i});
% %     pause;
% end
end

function []=PlotQuant(Fn,Nfn,leg)

Fn(Fn==0)=10^(-23);
Nfn(Nfn==0)=10^(-23);

[freq10 x10]=hist(log10(Fn),20);freq10=freq10./sum(freq10);
semilogx(10.^x10,freq10);
% mn_Fn=10^mean(log10(Fn));mx_Fn=max(freq10);
mn_Fn=10^mean(log10(Fn));mx_Fn=max(freq10);
hold on;
[freq10 x10]=hist(log10(Nfn),20);freq10=freq10./sum(freq10);
semilogx(10.^x10,freq10,'Color','r');
mn_Nfn=10^mean(log10(Nfn));mx_Nfn=max(freq10);
%%%% PLOT MEAN VALUES  ))))))))))))))))))
y_axis=(0:.1:1)*mx_Fn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Fn,y_axis,'b.-');
text(mn_Fn,mx_Fn,num2str(mn_Fn),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');
y_axis=(0:.1:1)*mx_Nfn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Nfn,y_axis,'r.-');
text(mn_Nfn,mx_Nfn,num2str(mn_Nfn),'BackgroundColor',[153 153 153]./256,'HorizontalAlignment','center','VerticalAlignment','baseline');

[p,h,s] = kruskalwallis([log10(Fn);log10(Nfn)],[ones(size(log10(Fn),1),1);ones(size(log10(Nfn),1),1)],'off');     display(p);
title_=['Kruskal Wallis p-value: ',num2str(p)];
title(title_);
ylabel('Percent of CEP Proteins (%)');
legend(leg,'Location','Best');

hold off;
end