function []=ExperimentalSProperties(IFiles,ExpIF)
close all;
[FileNames] = textread(IFiles,'%s',-1,'delimiter','\t');
file=length(FileNames);
Names=[];SELECTED=[];SPLEN=[];DESCR=[];TOTALL=[];MOLW=[];
for i=1:1:file
    [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(FileNames{i},'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
    Names=[Names;GnNames];
    SELECTED=[SELECTED;SEQ];
    SPLEN=[SPLEN;SPLen];
    DESCR=[DESCR;Descr];
    TOTALL=[TOTALL;TotLen];
    MOLW=[MOLW;MW];
end
[SubNames, Descr_sub, TotLen_sub, MW_sub, SPLen_sub, SEQ_sub] = textread(ExpIF,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

Peptides=length(SELECTED);
SelSeq_sp=cell(Peptides,1);
SelSeq_sp_sub=cell(length(SubNames),1);
for i=1:1:Peptides
    [AAs_int_sp SelSeq_sp(i)]=AARepresentation('int',1,SPLEN(i),{SELECTED{i}});
end
for i=1:1:length(SubNames)
    [AAs_int_sp_sub SelSeq_sp_sub(i)]=AARepresentation('int',1,SPLen_sub(i),{SEQ_sub{i}});
end

[Domains AllNames max_Nend max_Hend]=DomainsLoad(Names,'SignalP/SPDomains.txt');
[Domains_sub AllNames_sub max_Nend_sub max_Hend_sub]=DomainsLoad(SubNames,'SignalP/SPDomains.txt');

[posCharge Pho PhoSum pos_indx pho_indx]=CountPropertiesDomain(SelSeq_sp,AllNames,Names,Domains,max_Nend,max_Hend);

posCharge=posCharge(posCharge>=0);
Pho=Pho(Pho>=0);
PhoSum=PhoSum(PhoSum>=0);
Domains=Domains(posCharge>=0,:)

display('Experimental');
[posCharge_sub Pho_sub PhoSum_sub pos_indx_sub pho_indx_sub DomSeq_sub]=CountPropertiesDomain(SelSeq_sp_sub,AllNames_sub,SubNames,Domains_sub,max_Nend,max_Hend);

posCharge_sub=posCharge_sub(posCharge_sub>0);
Pho_sub=Pho_sub(Pho_sub>0);
PhoSum_sub=PhoSum_sub(PhoSum_sub>0);
Peptides=size(posCharge,1);
Peptides_sub=size(posCharge_sub,1);

index=(posCharge==0);
PrintSelection('Data/ZeroChargeNDomain.txt',{Names{index}},{DESCR{index}},TOTALL(index),MOLW(index),SPLEN(index),{SELECTED{index}});
index=(Domains(:,4)>=10);
PrintSelection('Data/LargeNDomain.txt',{Names{index}},{DESCR{index}},TOTALL(index),MOLW(index),SPLEN(index),{SELECTED{index}});

PrintProperties('Data/ExperimentalDomProp.txt',SubNames,DomSeq_sub,posCharge_sub,Pho_sub,Domains_sub(:,4),Domains_sub(:,5));

f_handle=figure(1);
x=1:1:max(max(Domains(:,4:end)));
[freq x]=hist(Domains(:,4),x);
freq=freq./Peptides;
    bar(x,freq,'y');hold on;
x=1:1:max(max(posCharge));
[freq x]=hist(posCharge,x);
freq=freq./Peptides;
    bar(x,freq,0.4,'r');hold off;
legend('N-domain Length',['Positive Charges (max: ',int2str(max(posCharge)),')'],'Location','best');
title(['Histogram of  N-Domain Lengths, Min-Max: ',int2str(min(Domains(:,4))),'-',int2str(max(Domains(:,4)))]);
xlabel('N-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle,'Figures/NDomainhist.bmp','bmp');

f_handle2=figure(2);
x=1:1:max(max(Domains(:,5:end)));
[freq x]=hist(Domains(:,5),x);
freq=freq./Peptides;
    bar(x,freq,0.8,'y');hold on;
x=1:1:max(max(Pho));
[freq x]=hist(Pho,x);
freq=freq./Peptides;
    bar(x,freq,0.4,'r');
[freq x]=hist(PhoSum);
freq=freq./Peptides;
    bar(x,freq,0.2,'g');hold off;
title(['Histogram of  H-Domain Lengths, Min-Max: ',int2str(min(Domains(:,5))),'-',int2str(max(Domains(:,5)))]);
xlabel('H-Domain length (aas)');
ylabel('Percentage of Proteins');
legend('H-domain Length',['Number of Hydrophobic Residues (max: ',int2str(max(Pho)),')'],'Engelman Hydrophobicity scale','Location','best');
saveas(f_handle2,'Figures/HDomainhist.bmp','bmp');


f_handle3=figure(3);
x=1:1:max(max(Domains(:,6:end)));
[freq x]=hist(Domains(:,6),x);
freq=freq./Peptides;
bar(x,freq,'y');
title(['Histogram of  C-Domain Lengths, Min-Max: ',int2str(min(Domains(:,6))),'-',int2str(max(Domains(:,6)))]);
xlabel('C-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle3,'Figures/CDomainhist.bmp','bmp');

f_handle4=figure(4);subplot(4,1,1);
x=1:1:max(max(posCharge));
[freq x]=hist(posCharge,x);
freq=freq./Peptides;
    bar(x,freq,'y');hold on;
[freq x]=hist(posCharge_sub,x);
freq=freq./Peptides_sub;
    bar(x,freq,0.4,'r');hold off;
legend('AllSecreted','Experimental Dataset','Location','best');
% title(['Histogram of  N-Domain Charges, Min-Max: ',int2str(min(posCharge)),'-',int2str(max(posCharge))]);
xlabel('N-Domain Charge (Instances of [KR])');
ylabel('Percentage of Proteins');
saveas(f_handle4,'Figures/NDomainCharge_compare.bmp','bmp');

f_handle5=figure(4);subplot(4,1,2);
x=1:1:max(max(Domains(:,4:end)));
[freq x]=hist(Domains(:,4),x);
freq=freq./Peptides;
    bar(x,freq,'y');hold on;
[freq x]=hist(Domains_sub(:,4),x);
freq=freq./Peptides_sub;
    bar(x,freq,0.4,'r');hold off;
% legend('AllSecreted','Experimental Dataset','Location','best');
% title(['Histogram of  N-Domain Length, Min-Max: ',int2str(min(Domains(:,4))),'-',int2str(max(Domains(:,4)))]);
xlabel('N-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle5,'Figures/NDomainLen_compare.bmp','bmp');

f_handle6=figure(4);subplot(4,1,3);
x=1:1:max(max(Pho));
[freq x]=hist(Pho,x);
freq=freq./Peptides;
    bar(x,freq,'y');hold on;
[freq x]=hist(Pho_sub,x);
freq=freq./Peptides_sub;
    bar(x,freq,0.4,'r');hold off;
% legend('AllSecreted','Experimental Dataset','Location','best');
% title(['Histogram of  H-Domain Hydrophobic Residues, Min-Max: ',int2str(min(Pho)),'-',int2str(max(Pho))]);
xlabel('H-Domain Hydrophobic Residues (Instances of [AGLMFIV])');
ylabel('Percentage of Proteins');
saveas(f_handle6,'Figures/HDomainPho_compare.bmp','bmp');

f_handle7=figure(4);subplot(4,1,4);
x=1:1:max(max(Domains(:,5:end)));
[freq x]=hist(Domains(:,5),x);
freq=freq./Peptides;
    bar(x,freq,'y');hold on;
[freq x]=hist(Domains_sub(:,5),x);
freq=freq./Peptides_sub;
    bar(x,freq,0.4,'r');hold off;
% legend('AllSecreted','Experimental Dataset','Location','best');
% title(['Histogram of  H-Domain Length, Min-Max: ',int2str(min(Domains(:,5))),'-',int2str(max(Domains(:,5)))]);
xlabel('H-Domain Length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle7,'Figures/HDomainLen_compare.bmp','bmp');
end

function []=PrintSelection(oFile,GnNames, Descr, TotLen, MW, SPLen, SEQ)

s=fopen(oFile,'w');
items=length(GnNames);
for i=1:1:items
    fprintf(s,'%s \t %s \t %d \t %d \t %d \t %s \t %s \t %s \n',GnNames{i},Descr{i},TotLen(i),MW(i),SPLen(i),SEQ{i},SEQ{i}(1:SPLen(i)),SEQ{i}(SPLen(i)+1:end));
end
fclose(s);

end

function []=PrintProperties(oFile,GnNames,DomSeq,posCharge,Pho,NDomain_len,HDomain_len)

s=fopen(oFile,'w');
items=length(posCharge);
fprintf(s,'Names \t N-Domain \t  H-Domain \t C-Domain \t Positively charged residues [KR] \t N-Domain Len (aas) \t H-Domain Hydrophobicity | H-Domain len \n');
for i=1:1:items
    fprintf(s,'%s \t %s \t %s \t %s \t %d \t %d \t %d | %d \n',GnNames{i},DomSeq{i,1},DomSeq{i,2},DomSeq{i,3},posCharge(i),NDomain_len(i),Pho(i),HDomain_len(i));
end
fclose all;

end
