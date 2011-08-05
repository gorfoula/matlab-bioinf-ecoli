function []=SignalPDomains(IFile,DFile)

[FileNames] = textread(IFile,'%s',-1,'delimiter','\t');
[DatasetFileNames] = textread(DFile,'%s',-1,'delimiter','\t');

posCharge=[];
Pho=[];
Domains=[];
Names={};
PhoSum=[];

for  i=1:1:length(FileNames)
    FileNames{i}
    [lines] = textread(FileNames{i},'%s',-1,'delimiter','\t');
    [Domains_cur Names_cur Catg]=ExtractDomainsSignalP(lines,int2str(i));
    [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(DatasetFileNames{i},'%s %s %d %d %d %s',-1,'delimiter','\t');  %%% read dataset file
    [posCharge_cur Pho_cur PhoSum_cur]=CountPropertiesDomain(SEQ,GnNames,Names_cur,Domains_cur);
    posCharge=[posCharge_cur;posCharge];
    Pho=[Pho_cur;Pho];
    Names_cur=[Names_cur Catg];
    Domains=[Domains_cur;Domains];
    Names=[Names_cur;Names];
    PhoSum=[PhoSum_cur;PhoSum];
end

largeND_nms=Names(Domains(:,4)>=13,:);
largeHD_nms=Names(Domains(:,5)>=13,:);
%%%%%%%% FIGURES AND SAVE  %%%%%%%%%%
s=fopen(['SignalP/','LargeNDomains.txt'],'w');
s2=fopen(['SignalP/','LargeHDomains.txt'],'w');
s3=fopen(['SignalP/','SPDomains.txt'],'w');

for i=1:1:length(largeND_nms(:,1))
    switch largeND_nms{i,2}
        case '1'
            catg='Periplasmic';
        case '2'
            catg='outer Membrane b-barrel';
        case '3'
            catg='outer Membrane lipoproteins';
    end
    fprintf(s,'%s\t%s\n',largeND_nms{i,1},catg);
end
for i=1:1:length(largeHD_nms(:,1))
    switch largeHD_nms{i,2}
        case '1'
            catg='Periplasmic';
        case '2'
            catg='outer Membrane b-barrel';
        case '3'
            catg='outer Membrane lipoproteins';
    end
    fprintf(s2,'%s\t%s\n',largeHD_nms{i,1},catg);
end
% fprintf(s3,'Gene Name\tN-Domain end position\tH-Domain end position\tC-Domain end position\tN-Domain length\tH-Domain length\tC-Domain length\n');
for i=1:1:length(Domains(:,1))
    fprintf(s3,'%s\t',Names{i,1});
    fprintf(s3,'%d\t%d\t%d\t%d\t%d\t%d\n',Domains(i,1),Domains(i,2),Domains(i,3),Domains(i,4),Domains(i,5),Domains(i,6));
end
fclose all;
close all;

Peptides=length(Domains(:,1));

f_handle=figure(1);
x=1:1:max(max(Domains(:,4:end)));
[freq x]=hist(Domains(:,4),x);
freq=freq./Peptides;
bar(x,freq);
x=1:1:max(max(posCharge));
[freq x]=hist(posCharge,x);
freq=freq./Peptides;
hold on;bar(x,freq,0.4,'r');
legend('N-domain Length',['Positive Charges (max: ',int2str(max(posCharge)),')'],'Location','best');
title(['Histogram of  N-Domain Lengths, Min-Max: ',int2str(min(Domains(:,4))),'-',int2str(max(Domains(:,4)))]);
xlabel('N-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle,'Figures/NDomainhist.bmp','bmp');

f_handle2=figure(2);
x=1:1:max(max(Domains(:,5:end)));
[freq x]=hist(Domains(:,5),x);
freq=freq./Peptides;
bar(x,freq,0.8);
x=1:1:max(max(Pho));
[freq x]=hist(Pho,x);
freq=freq./Peptides;
hold on;bar(x,freq,0.4,'r');
[freq x]=hist(PhoSum);
freq=freq./Peptides;
hold on;bar(x,freq,0.2,'g');
legend('H-domain Length',['Number of Hydrophobic Residues (max: ',int2str(max(Pho)),')'],'Engelman Hydrophobicity scale','Location','best');
title(['Histogram of  H-Domain Lengths, Min-Max: ',int2str(min(Domains(:,5))),'-',int2str(max(Domains(:,5)))]);
xlabel('H-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle2,'Figures/HDomainhist.bmp','bmp');

f_handle3=figure(3);
x=1:1:max(max(Domains(:,6:end)));
[freq x]=hist(Domains(:,6),x);
freq=freq./Peptides;
bar(x,freq);
title(['Histogram of  C-Domain Lengths, Min-Max: ',int2str(min(Domains(:,6))),'-',int2str(max(Domains(:,6)))]);
xlabel('C-Domain length (aas)');
ylabel('Percentage of Proteins');
saveas(f_handle3,'Figures/CDomainhist.bmp','bmp');

end