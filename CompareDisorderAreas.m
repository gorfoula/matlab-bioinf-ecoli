function []=CompareDisorderAreas(DatasetFiles,DisorderFiles)

close all;

[DataFileNames] = textread(DatasetFiles,'%s',-1,'delimiter','\t');
[DisProtFileNames] = textread(DisorderFiles,'%s',-1,'delimiter','\t');

[disorderRegions_cyto count_cyto total_cyto]=ReadDisorderAreas(DataFileNames{1},DisProtFileNames{1});
[disorderRegions count total]=ReadDisorderAreas(DataFileNames{2},DisProtFileNames{2});

found=regexp(DisProtFileNames{1},'[_/]');

len=[min([disorderRegions(:,1);disorderRegions_cyto(:,1)]) max([disorderRegions(:,1);disorderRegions_cyto(:,1)])];
position=[min([disorderRegions(:,2);disorderRegions_cyto(:,2)]) max([disorderRegions(:,2);disorderRegions_cyto(:,2)])];
total_regns=[min([total';total_cyto']) max([total';total_cyto'])];

f1=figure(1);
[freq x]=hist(disorderRegions(:,2),5:10:position(2));
bar(x,freq./count,'r');
hold on;
[freq x]=hist(disorderRegions_cyto(:,2),5:10:position(2));
bar(x,freq./count_cyto,0.4,'c');
xlabel('Position on Mature Sequence');
ylabel('Precentage over all disorder regions');
title(['Histogram of the position disordered regions on Mature domain (DisProt ',DisProtFileNames{1}(found(1)+1:found(2)-1),')']);
legend('Secreted','Cytoplasmic','Location','best');
f2=figure(2);
[freq x]=hist(disorderRegions(:,1),5:10:len(2));
bar(x,freq./count,'r');
hold on;
[freq x]=hist(disorderRegions_cyto(:,1),5:10:position(2));
bar(x,freq./count_cyto,0.4,'c');
xlabel('Length of disorder region');
ylabel('Precentage over all disorder regions');
title(['Histogram of the length of disordered regions on Mature domain (DisProt ',DisProtFileNames{1}(found(1)+1:found(2)-1),')']);
legend('Secreted','Cytoplasmic','Location','best');
f3=figure(3);
[freq x]=hist(disorderRegions(:,1),total_regns(1):1:total_regns(2));
bar(x,freq./count,'r');
hold on;
[freq x]=hist(total_cyto,total_regns(1):1:total_regns(2));
bar(x,freq./count_cyto,0.4,'c');
xlabel('Number of disorder regions per protein');
ylabel('Precentage over all proteins');
title(['Histogram of the number of disordered regions on Mature domain (DisProt ',DisProtFileNames{1}(found(1)+1:found(2)-1),')']);
legend('Secreted','Cytoplasmic','Location','best');


saveas(f1,['Figures/DisorderRegPosHist_',DisProtFileNames{1}(found(1)+1:found(2)-1),'.bmp'],'bmp');
saveas(f2,['Figures/DisorderRegLenHist_',DisProtFileNames{1}(found(1)+1:found(2)-1),'.bmp'],'bmp');
saveas(f3,['Figures/DisorderNumRegHist_',DisProtFileNames{1}(found(1)+1:found(2)-1),'.bmp'],'bmp');


end