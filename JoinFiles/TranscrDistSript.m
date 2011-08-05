function []=TranscrDistSript(if_PreomicsLists,if_ProteinLevel,if_mRNALevel,sel)
close all;
% Merge_mRNA_moleculesPerCell.txt
%Merge_Protein_moleculesPerCell.txt

[TABLE index_wide files_wide print_cols_wide]=ExludeCyto(if_PreomicsLists,2,0);
[Check_Table]=CellTable2Double(CELLENV_wide_studies(2:end,files_wide_studies+2));
NotDetectedProtLevel=(Check_Table==0);

str_point=6;

[CELLENV_mRNA NotDetMRNA Detected_Us EmptyCells_mRNA LowDet MBE]=PlotFounVsNotFound(if_mRNALevel,2,sel,str_point);
[CELLENV_P NotDetProt Detected_Us EmptyCells_P LowDet_P MBE_P USDETECTED USNOTDETECTED Quant]=PlotFounVsNotFound(if_ProteinLevel,2,sel,str_point);

detprot_not_mRNA=((NotDetMRNA) & not(EmptyCells_mRNA) & not(NotDetProt));
figure;PlotQuant(Quant(USDETECTED(:,7)& detprot_not_mRNA,7),Quant(USNOTDETECTED(:,7) & detprot_not_mRNA,7),[{'All CEP'};{'All CEP'};{'NOT Found by us'};{'NOT Found by us'}]);

Abudance=cell(4407,1);
Abudance(NotDetMRNA,:)={'A'};
Abudance(LowDet,:)={'L'};
Abudance(EmptyCells_mRNA,:)={'E'};
Abudance(not(EmptyCells_mRNA | LowDet |NotDetMRNA),:)={'P'};

FileWriteTable('JoinFiles/Transcriptomics/AllmRNAinfo.txt',[CELLENV_mRNA [{'Abundance Level'};Abudance]],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/mRNAEmpty.txt',CELLENV_mRNA(logical([1;and(EmptyCells_mRNA,Detected_Us)]),:),[],'w');
% FileWriteTable('JoinFiles/Transcriptomics/absent.txt',CELLENV_mRNA(logical([1;and(NotDetMRNA,Detected_Us)]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/absent.txt',CELLENV_mRNA(logical([1;NotDetMRNA]),:),[],'w');

FileWriteTable('JoinFiles/Transcriptomics/mRNA_present.txt',CELLENV_mRNA(logical([1;and(and(not(NotDetMRNA),not(EmptyCells_mRNA)),Detected_Us)]),:),[],'w');

FileWriteTable('JoinFiles/Transcriptomics/All_Protein_info.txt',CELLENV_P,[],'w');
FileWriteTable('JoinFiles/Transcriptomics/Protein_Empty.txt',CELLENV_P(logical([1;and(not(EmptyCells_P),Detected_Us)]),:),[],'w');
FileWriteTable('JoinFiles/Transcriptomics/Protein_absent.txt',CELLENV_P(logical([1;and(and(not(NotDetProt),not(EmptyCells_P)),Detected_Us)]),:),[],'w');

% selcol=8;
% MolPerCell=CellTable2Double(CELLENV_P(2:end,selcol));MolPerCell=MolPerCell(not(EmptyCells_P));
% MolPerCell_log=log10(MolPerCell+1);
% BINS=100;
% g=min(MolPerCell_log):(max(MolPerCell_log)-min(MolPerCell_log))/(BINS-1):max(MolPerCell_log);
% g=(10.^g)-1;
% for j=1:1:BINS 
%     DESCR_MolPerCell=zeros(length(MolPerCell),1);
%     DESCR_MolPerCell(MolPerCell>g(j))=1;
%     [table,X2(j),pvalue_X2(j)] = crosstab(Check_Table(not(EmptyCells_P)), DESCR_MolPerCell);
% %     [pvalue_X2(j) X2(j)] = independenceTest(Check_Table, DESCR_MolPerCell, [max(Check_Table) 2], 1);
% end
% 
% [freq x]=hist(MolPerCell,g);
% figure;
% semilogx(x,freq*100./sum(freq),'LineWidth',2);hold on;
% semilogx(x,pvalue_X2,'Color','r','LineWidth',2);
% legend([{CELLENV_P{1,selcol}};{'p-value independence'}]);
% xlabel(CELLENV_P{1,selcol});
% ylabel({'P-value' 'Percent of Proteins'});
% hold off;
% figure;
% bar(g,freq*100./sum(freq),'LineWidth',2);hold on;
% plot(x,pvalue_X2,'Color','r','LineWidth',2);
% xlabel(CELLENV_P{1,selcol});
% ylabel({'P-value' 'Percent of Proteins'});
% hold off;


Empty=EmptyCells_P & EmptyCells_mRNA & NotDetectedProtLevel;

%[pvalue_X2 X2] = independenceTest(Check_Table, DESCR_MolPerCell' , [max(Check_Table) 2], 0);

Not_detected_at_any_level=NotDetProt & NotDetMRNA & not(Detected_Us);
Not_detected_at_mRNA_or_P=(NotDetProt | NotDetMRNA) & not(Detected_Us);
NotMRNA_but_Prot=not(NotDetProt) & NotDetMRNA & not(Detected_Us);
Detected_at_least_one_level=(not(NotDetProt) | not(NotDetMRNA)) & not(Detected_Us);
Detected_at_transcr=not(NotDetMRNA) & not(Detected_Us);
Detected_at_protein=not(NotDetProt) & not(Detected_Us);

display(['Empty: ',num2str(sum(Empty))]);
display(['NoTranscript: ',num2str(sum(NotDetMRNA & not(Detected_Us)))]);
display(['NoProteinLevel: ',num2str(sum(NotDetProt & not(Detected_Us)))]);
display(['Not_detected_at_any_level: ',num2str(sum(Not_detected_at_any_level))]);
display(['Not_detected_at_mRNA_or_P: ',num2str(sum(Not_detected_at_mRNA_or_P))]);
display(['NotMRNA_but_Prot: ',num2str(sum(NotMRNA_but_Prot))]);
display(['Detected_at_transcript: ',num2str(sum(Detected_at_transcr))]);
display(['Detected_at_protein: ',num2str(sum(Detected_at_protein))]);
display(['Detected_at_least_one_level: ',num2str(sum(Detected_at_least_one_level))]);

FileWriteTable('JoinFiles/Transcriptomics/NotDetectedAtAnyLevel.txt',[CELLENV_mRNA(logical([1;Not_detected_at_any_level]),:) CELLENV_P(logical([1;Not_detected_at_any_level]),2:end)],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/NotDetectedOnlyTransr.txt',[CELLENV_mRNA(logical([1;NotMRNA_but_Prot]),:) CELLENV_P(logical([1;NotMRNA_but_Prot]),2:end)],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/DetectedBothLevels.txt',[CELLENV_mRNA(logical([1;Detected_at_least_one_level]),:) CELLENV_P(logical([1;Detected_at_least_one_level]),2:end)],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/DetectedProtein.txt',[CELLENV_mRNA(logical([1;Detected_at_protein]),:) CELLENV_P(logical([1;Detected_at_protein]),2:end)],[],'w');
FileWriteTable('JoinFiles/Transcriptomics/DetectedTranscr.txt',[CELLENV_mRNA(logical([1;Detected_at_transcr]),:) CELLENV_P(logical([1;Detected_at_transcr]),2:end)],[],'w');

FileWriteTable('JoinFiles/Transcriptomics/Empty.txt',[CELLENV_mRNA(logical([1;Empty]),:) CELLENV_P(logical([1;Empty]),2:end)],[],'w');

% subplot(2,1,1);PlotQuant(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),leg);xlabel(Header{i});
% subplot(2,1,2);SimpleHist(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),[leg(1) leg(3)]);xlabel(Header{i});

end

function [CELLENV NotDet Detected_Us EmptyCells LowDet MBE USDETECTED USNOTDETECTED Quant]=PlotFounVsNotFound(FilePath,ref_file,sel,str_point)

%%% sel:    1 for comparing found not found
%%%         2 for cmparing all protein with all in last file list
%%%         3 rerun (do not read previous saved CEP files) compare found VS not found
%%%         4 rerun (do not read previous saved CEP files)  and compare All proteins with all in last file list

if(sel>2)
    [CELLENV index files]=ExludeCyto(FilePath,ref_file,1);
else
    [CELLENV index files]=ExludeCyto(FilePath,ref_file,0);
end

% CELLENV=[CELLENV;CYTO];

[MobileElem]=CellTable2Double(CELLENV(2:end,files+str_point-1));

[Quant index]=CellTable2Double(CELLENV(2:end,files+str_point:end-1));
[Check_Table]=CellTable2Double(CELLENV(2:end,2:files+2));
Header=CELLENV(1,files+str_point:end-1);
Detected_Us=Check_Table(:,files)==1;
%  (DETECTED & (Measured At Zero | Not Measured)) | NOT DETECTED
EmptyCells=(sum(index,2)==0) | Check_Table(:,end)==0;
NotDet=( (Check_Table(:,end)>0) & (sum(Quant(:,2:end-1),2)<=1) );
LowDet=( NotDet & ((Quant(:,end)<=1) & index(:,end)) );
MBE=(MobileElem==1) & (NotDet==1);

% Detected_Us_square = ones(size(Detected_Us'))'*Detected_Us';
Detected_Us_square = Detected_Us*ones(size(Detected_Us'));
switch(sel)
    case {1,3}
        leg=[{'Found by us'};{'Found by us'};{'NOT Found by us'};{'NOT Found by us'}];
        USDETECTED=logical(Detected_Us_square(:,1:size(index,2)).*index);
        USNOTDETECTED=logical(not(Detected_Us_square(:,1:size(index,2))).*index);
    case {2,4}
        leg=[{'All CEP'};{'All CEP'};{'NOT Found by us'};{'NOT Found by us'}];
        USDETECTED=index;
        USNOTDETECTED=logical(Detected_Us_square(:,1:size(index,2)).*index);
    otherwise
        display('Sel variable should be 1 or 2');
end

% Quant=Quant+1;

for i=1:length(Header)
    f=figure;
    subplot(2,1,1);PlotQuant(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),leg);xlabel(Header{i});
    subplot(2,1,2);SimpleHist(Quant(USDETECTED(:,i),i),Quant(USNOTDETECTED(:,i),i),[leg(1) leg(3)]);xlabel(Header{i});
end
f=figure;

for k=1:1:size(Quant,1)
    A=zeros(size(Quant,1),1);
    B=zeros(size(Quant,1),1);
    if(sum(USDETECTED(k,:))>0)
        A(k)=round( sum(Quant(k,USDETECTED(k,:)))./sum(USDETECTED(k,:)) );
    end
    if(sum(USNOTDETECTED(k,:))>0)
        B(k)=round( sum(Quant(k,USNOTDETECTED(k,:)))./sum(USNOTDETECTED(k,:)) );
    end
end
SimpleHist(A,B,[leg(1) leg(3)]);xlabel('Normalized mRNA Level');

end

function []=PlotQuant(Fn,Nfn,leg)

[freq10 x10]=hist(log10(Fn),20);freq10=freq10./sum(freq10);
semilogx(10.^x10,freq10);
mn_Fn=10^mean(log10(Fn));mx_Fn=max(freq10);
hold on;
bar(10.^x10,freq10);
[freq10 x10]=hist(log10(Nfn),20);freq10=freq10./sum(freq10);
semilogx(10.^x10,freq10,'Color','r');
mn_Nfn=10^mean(log10(Nfn));mx_Nfn=max(freq10);
h=bar(10.^x10,freq10,0.7,'r');

ylabel('Percent of CEP Proteins (%)');
legend(leg);

y_axis=(0:.1:1)*mx_Fn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Fn,y_axis,'b.');
text(mn_Fn,mx_Fn,num2str(mn_Fn));
y_axis=(0:.1:1)*mx_Nfn;
x_axis=ones(length(y_axis),1);
plot(x_axis*mn_Nfn,y_axis,'r.');
text(mn_Nfn,mx_Nfn,num2str(mn_Nfn));

[p,h,s] = kruskalwallis([log10(Fn);log10(Nfn)],[ones(size(log10(Fn),1),1);ones(size(log10(Nfn),1),1)],'off');     display(p);
title_=['Kruskal Wallis p-value: ',num2str(p)];
title(title_);

hold off;

end

function []=SimpleHist(Fn,Nfn,leg)

if(length(unique(Fn))==1)
    [f_f x_f]=hist(Fn);f_f=f_f./sum(f_f);
    [f_n x]=hist(Nfn);f_n=f_n./sum(f_n);
else
    [f_f x_f]=hist(Fn,unique(Fn));f_f=f_f./sum(f_f);
    [f_n x]=hist(Nfn,x_f);f_n=f_n./sum(f_n);
end
bar(x_f,f_f);
hold on;
bar(x_f,f_n,0.7,'r');
hold off;
ylabel('Percent of CEP Proteins (%)');
legend(leg);

[p,h,s] = kruskalwallis([Fn;Nfn],[ones(size(Fn,1),1);ones(size(Nfn,1),1)],'off');     display(p);
title_=['Kruskal Wallis p-value: ',num2str(p)];
title(title_);

end