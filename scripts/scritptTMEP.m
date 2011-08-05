
Proteins=[13	39	886	161	1140	358	1469	1366	574	56	54	2049	295	104	631	299	6	73	181	435	42	22	213	81	1659	1802	64	42	2599	405	54	72	633 132];
Names=[{'Molloy1999'} {'Molloy2000'} {'Gevaert2002'} {'Yan2002'} {'Corbin2003'} {'Fountoulakis2003'} {'Taoka2004'} {'Butland2005'} {'Lopez_Campistrous2005'} {'Spelbrink2005'} {'Stenberg2005'} {'Arifuzzaman2006'} {'Baars2006'} {'Huang2006'} {'Ji2006'} {'Lasserre2006'} {'Marani2006'} {'Cirulli2007'} {'Wagner2007'} {'Zhang2007'} {'Jarchow2008'} {'Qian2008'} {'Wagner2008'} {'Xia2008'} {'Iwasaki2009'} {'Masuda2009'} {'Vertommen2009'} {'Hemm2010'} {'Iwasaki2010'} {'Muller2010'} {'Pan2010'} {'Price2010'} {'Thein'} {'Wickstrom2010'}];

x=0:100:max(Proteins);
figure;hist(Proteins,x);
ylabel('#studies');
xlabel('#Proteins');


TheoReadTable_=ReadTable('JoinFiles/Transcriptomics/Masuda2009_fin2.txt');
TheoReadTable_=TheoReadTable_(2:end,:);

[Quant index]=CellTable2Double(TheoReadTable_(:,8));

log_10=log10(Quant+1);
figure;[freq x]=hist(log_10(index));bar(x,freq);
mx=x(find(freq==max(freq)));
mn=mean(log_10(index));
sigma=std(log_10(index));
th1=mx-(sigma/2);th2=mx+(sigma/2);
low=(log_10<=th1) & index;
hight=(log_10>=th2) & index;
mid=not(low) & not(hight);

Levels=zeros(length(index),1);
Levels(hight)=2;
Levels(mid)=1;
% 
% TheoIF='JoinFiles/TheoreticalPeptides/NotDetected.txt';
% PepTheo_if='JoinFiles/TheoreticalPeptides/NotDetected_PepTheo.txt';
% CountTheo_if='JoinFiles/TheoreticalPeptides/NotDetected_IDs_Theo.txt';
TheoIF='JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011CEPNonDet.txt';
PepTheo_if='JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011CEPNonDet_PepTheo.txt';
CountTheo_if='JoinFiles/OtherStudies/Notepad Lists/OurStudy_Feb2011CEPNonDet_IDs_Theo.txt';
col=[1 6];

[UniqueIDs_Theo UNIQUE_PEP_Theo]=CountLineSameID(TheoIF,[col(1) col(2)],1);
FileWriteTable(PepTheo_if,UNIQUE_PEP_Theo,[],'w');
FileWriteTable(CountTheo_if,UniqueIDs_Theo,[],'w');