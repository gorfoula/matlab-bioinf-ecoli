function []=ExpResAnalysis(DatasetF,DomainsF,ExpResF)

[GnNames_exp, Descr, TotLen, MW, SPLen, SEQ] = textread(DatasetF,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[AllNames, ND_end, HD_end, CD_end, ND_len, HD_len, CD_len] = textread(DomainsF,'%s %d %d %d %d %d %d',-1,'delimiter','\t');  %%% read dataset file
[Names SecrEff Refolding] = textread(ExpResF,'%s %f %f',-1,'delimiter','\t');  %%% read dataset file


[positives Pho PhoSum]=CountPropertiesDomain(SEQ,AllNames,GnNames_exp,[ND_end, HD_end, CD_end, ND_len, HD_len, CD_len])



index=(and(Refolding<1.5,not(Refolding==0)))

table=[SecrEff(index) Refolding(index) positives(index) Pho(index) PhoSum(index)];

[RHO PV]=corr(table,'type','Pearson');
RHO=tril(RHO);

features={'Secretion Efficiency' 'Refolding' 'Positive Charged aas' 'Hydrophobic aas' 'Hydrophobi scale'};
%%%%%%%%%%%%
index=(and(Refolding>=1.5,not(Refolding==0)))

table=[SecrEff(index) Refolding(index) positives(index) Pho(index) PhoSum(index)];

[RHO_h PV_h]=corr(table,'type','Pearson');
RHO_h=tril(RHO_h);

%%%%%%%%%%%%

s=fopen('Data/LowRefolding.txt','w');
s2=fopen('Data/HighRefolding.txt','w');
fprintf(s,'\tSecretion Efficiency\tRefolding\tPositive Charged aas\tHydrophobic aas\tHydrophobi scale\n');
fprintf(s2,'\tSecretion Efficiency\tRefolding\tPositive Charged aas\tHydrophobic aas\tHydrophobi scale\n');
for i=1:1:length(features)
    fprintf(s,'%s\t',features{i});
    fprintf(s2,'%s\t',features{i});
    for j=1:1:i-1
        fprintf(s,'%1.3f|%1.3f\t',RHO(i,j),PV(i,j));
        fprintf(s2,'%1.3f|%1.3f\t',RHO_h(i,j),PV_h(i,j));
    end
    fprintf(s,'\n');
    fprintf(s2,'\n');
end

fclose all;

hist(Refolding(Refolding>0));
end