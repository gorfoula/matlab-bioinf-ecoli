close all;clear all; 

% file='JoinFiles/Peripheral/List_potential_F1.txt';
% [FOUNDTABLE]=FishLines(file,'JoinFiles/Peripheral/ECOLI_K12_4407_TMHMM_v3.txt',[1 0],[1 3],[1],1);
% [FOUNDTABLE]=FishLines(file,'JoinFiles/Peripheral/ECOLI_K12_4407_PHOBIOUS_v3.txt',[1 0],[1 3],[],1);
% file1='JoinFiles/Peripheral/InnerMem_NoTMs_Potential_F1.txt'
% [FOUNDTABLE]=FishLines(file1,'JoinFiles/Peripheral/ECOLI_K12_4407_TMHMM_v3.txt',[1 0],[1 3],[1],1);
% [FOUNDTABLE]=FishLines(file1,'JoinFiles/Peripheral/ECOLI_K12_4407_PHOBIOUS_v3.txt',[1 0],[1 3],[],1);



[CATCH_] = ReadTable('JoinFiles/Peripheral/ECOLI_K12_4407_TMHMM_v3.txt','\n');

CATCH_=CATCH_(:,3:end);

% LowTms=CellTable2Double(CATCH_(:,3))<.9;
CATCH=CATCH_(CellTable2Double(CATCH_(:,3))<1,:);


[TABLE Address]=Douplicates('JoinFiles/Peripheral/ECOLI_K12_4407_TMHMM_v3.txt',[0 5],CATCH);



% [CATCH_all] = ReadTable('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','\n');
% [TABLE_all Address_all]=Douplicates('JoinFiles/k12/TableS3_K12_ProteomeV13.txt',[1 5],CATCH_all);

% 
% CATCH=CATCH(LowTMs,:);
% [TABLE Address]=Douplicates('JoinFiles/Peripheral/TMHMM_mgto_AF1B.txt',[0 4],CATCH);

all_B=(Address==find(strcmp(TABLE(:,1),'B')));
all_F1=(Address==find(strcmp(TABLE(:,1),'F1')));
all_A=(Address==find(strcmp(TABLE(:,1),'A')));
all_PF1=(Address==find(strcmp(TABLE(:,1),'Potential F1')));
all_PF1B=(Address==find(strcmp(TABLE(:,1),'Potential F1 (B)')));

withTMs=CellTable2Double(CATCH(:,2))>0;

BwithTMareas=withTMs & all_B;
F1withTMareas=withTMs & all_F1;
AwithTMareas=withTMs & all_A;
PF1withTMareas=withTMs & all_PF1;
PF1BwithTMareas=withTMs & all_PF1B;

[TABLE_B_all Address_B_all]=Douplicates('..txt',[0 1],CATCH(all_B,:));
[TABLE_F1_all Address_F1_all]=Douplicates('..txt',[0 1],CATCH(all_F1,:));
[TABLE_A_all Address_A_all]=Douplicates('..txt',[0 1],CATCH(all_A,:));
[TABLE_PF1_all Address_PF1_all]=Douplicates('..txt',[0 1],CATCH(all_PF1,:));
[TABLE_PF1B_all Address_PF1B_all]=Douplicates('..txt',[0 1],CATCH(all_PF1B,:));

[TABLE_B Address_B]=Douplicates('..txt',[0 1],CATCH(BwithTMareas,:));
[TABLE_F1 Address_F1]=Douplicates('..txt',[0 1],CATCH(F1withTMareas,:));
[TABLE_A Address_A]=Douplicates('..txt',[0 1],CATCH(AwithTMareas,:));
[TABLE_PF1 Address_PF1]=Douplicates('..txt',[0 1],CATCH(PF1withTMareas,:));
[TABLE_PF1B Address_PF1B]=Douplicates('..txt',[0 1],CATCH(PF1BwithTMareas,:));

% CellTable2Double
perc_ov_all=Double2CellTable(ceil([max(Address_B)/max(Address_B_all);max(Address_F1)/max(Address_F1_all);max(Address_A)/max(Address_A_all);max(Address_PF1)/max(Address_PF1_all);max(Address_PF1B)/max(Address_PF1B_all)]*100),1);

% [Names Address_of_Names]=Douplicates('..txt',[0 1],CATCH(:,:));
[K12] = ReadTable('JoinFiles/k12/TableS3_K12_ProteomeV13.txt','\n');
[FOUNDTABLE]=FishLines('.txt','2.txt',[1 0],[1 1],[11],1,K12,CATCH(withTMs,1));

groups=Address(withTMs)-1;
categoty=MergeColumns([TABLE(:,1) perc_ov_all],' %over all: ');
temp=CATCH(withTMs,2);
Y=CellTable2Double(temp);
x=min(Y):5:50;
FigureLegends(x,Y,groups,[{'Predicted TMs'} {'Len'} {'Percent over all areas'}],categoty,'p',{'-','.';':','o';'-','.';'-','.';'-','.'},0,[1 3 1 1]);

Y_pos=CellTable2Double(CATCH(withTMs,4))./CellTable2Double(FOUNDTABLE(2:end,1));
Y2=CellTable2Double(CATCH(withTMs,3));
Y_log=log10(Y2);
x_=min(Y_log):(max(Y_log)-min(Y_log))/150:max(Y_log);
x=0:.01:1;
[Colors]=FigureLegends(x,Y2,groups,[{' '} {'Max Prob'} {'Percent over all areas'}],categoty,'p',{'-','.';':','o';'-','.';'-','.';'-','.'},0,[1 3 1 2]);
% figure;semilogx(10.^x_,hist(Y(groups==0),10.^x_),'Color',Colors(1,:));hold on;
% semilogx(10.^x_,hist(Y(groups==1),10.^x_),'Color',Colors(2,:));
% semilogx(10.^x_,hist(Y(groups==2),10.^x_),'Color',Colors(3,:));hold off;
% legend(categoty);

Y3=CellTable2Double([TABLE_B(:,2);TABLE_F1(:,2);TABLE_A(:,2);TABLE_PF1(:,2)]);
groups2=[zeros(size(TABLE_B,1),1);ones(size(TABLE_F1,1),1)*1;ones(size(TABLE_A,1),1)*2;ones(size(TABLE_PF1,1),1)*3];
x=min(Y3):1:max(Y3);
FigureLegends(x,Y3,groups2,[{' '} {'#Areas per protein'} {'Percent over all proteins with such areas'}],categoty,'p',{'-','.';':','o';'-','.';'-','.';'-','.'},0,[1 3 1 3]);

% close all;
indx=(groups<3);
r=.008;
figure;hold on;
xlim([0 1]);ylim([0 1]);
Myscatter(Y_pos(indx)',Y2(indx)',groups(indx),r);
legend('Cytoplasmic','Inner Membrane','Peripheral');
ylabel('Max Probability');
xlabel('Pos');hold off;

figure;hold on;
Myscatter(Y(indx)',Y2(indx)',groups(indx),r);
legend('Cytoplasmic','Inner Membrane','Peripheral');
ylabel('Max Probability');
xlabel('TM lenght');hold off;
xlim([0 1]);
ylim([0 max(Y(indx))]);








% h1=scatter(Y2(groups==0),Y_pos(groups==0));
% h2=scatter(Y2(groups==1),Y_pos(groups==1));
% h3=scatter(Y2(groups==2),Y_pos(groups==2));
% set(h3,'Marker','.','SizeData',72);
% set(h1,'Marker','.','SizeData',72);
