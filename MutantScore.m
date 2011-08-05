function []=MutantScore(MutantFile,SelFeatFILE,WeigthFile,len,sel,VALUES,nplot)

[mutantfile]=SynthesizeMutantFile(MutantFile,'Dataset/Allsecreted.txt');

[Name, SP, Mature, FolwSeq, Efficiency, Ref, Comments] = textread(MutantFile,'%s %s %s %s %d %s %s',-1,'delimiter','\t');
[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(mutantfile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[aa, FeatureID] = textread(SelFeatFILE,'%d %d',-1,'delimiter','\t');
FeatureID=FeatureID-1;
% VALUES=20;

if(exist('len','var')==0)
    NumericCode=mod(max(FeatureID),VALUES);
    PosOnPepseq=floor(max(FeatureID)/VALUES)+(NumericCode>0);
    len=PosOnPepseq;
end

Efficiency=Efficiency(TotLen>=SPLen+len);
MOLECULARW=MW(TotLen>=SPLen+len);
DESCRIPTION=Descr(TotLen>=SPLen+len);
Names=GnNames(TotLen>=SPLen+len);   %Gene Names
SELECTED=SEQ(TotLen>=SPLen+len);    %Select Petides with Seq long enough
samples=length(SELECTED);
TOTALLEN=TotLen(TotLen>=SPLen+len);
SPLen=SPLen(TotLen>=SPLen+len);     %Select corresponding Signal Peptide Lengths

WT_SP_score=zeros(samples,1);
AAs3D=zeros(samples,VALUES,len);
for peptide=1:1:samples
    [AAs3D(peptide,:,:) SelSeq(peptide)]=AARepresentation('bin',2,len+1,{SELECTED{peptide}});  %% Methionine excluded
end  
if(sel==20119)
    [CatgNames,Catg3D]=CatgRepr(AAs3D,'bin','Features/mature_AAProperties.txt');
    [CatgNames_9f,Catg3D_9f]=CatgRepr(AAs3D,'bin','Features/AAProperties_new.txt');
    Catg3D=reshape(Catg3D,samples,len*11);
    Catg3D_9f=reshape(Catg3D_9f,samples,len*9);
    AAs3D=reshape(AAs3D,samples,VALUES*len);  
    AAs3D=[AAs3D Catg3D Catg3D_9f];
else
    AAs3D=reshape(AAs3D,samples,VALUES*len);  
end
    
    [w] = textread(WeigthFile,'%f',-1,'delimiter','\t');
    score=AAs3D(:,FeatureID)*w(1:end-1)+w(end);
    
    groups=Efficiency>0;
    cp = classperf(double(groups));
    classes=score>0;
    classperf(cp,double(classes));
    cp.CorrectRate
    cp.ErrorRate
    
    score=score./(score(1))*100;
    score(score<0)=0;
    PlotMutantScores(Names,score,samples,Efficiency,nplot);
end

function []=PlotMutantScores(Names,score,samples,Efficiency,nplot)
if(nplot==1)
    close all;
else
    hold on;
end
Red=rand(54,1);Green=0:.01:1;Blue=rand(length(Red),1);
f1=figure(1);subplot(2,1,nplot);
x=1:1:samples;
bar(x,score,'y');hold on;
bar(x,Efficiency,0.4,'r');
ylim([0 max([Efficiency;score])+10]);
count=1;
for j=1:1:samples
    even=mod(count,2);
    if(even)
        text(x(j),Efficiency(j)+1,Names{j},'Color','k','FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');
    else
        text(x(j),Efficiency(j)-1,Names{j},'Color','k','FontUnits','points','FontSize',10,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','top');
    end
    count=count+1;
end
legend('Model','Experimental','Location','best');
title('Experimental Secretion Measurements of phoAKL Mutants VS Model estimation');
xlabel('Mutant a/a');
ylabel('Secretion/score (persecent of phoA secretion)');
grid on;
hold off;
end