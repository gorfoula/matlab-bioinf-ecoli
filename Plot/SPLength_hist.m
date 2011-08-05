function []=SPLength_hist(PeriIF,bBarrelIF,LipoIF)

[P_names, P_Descr, P_TotLen, P_MW, P_SPLen, P_TOTseq] = textread(PeriIF,'%s %s %d %d %d %s',-1,'delimiter','\t');
[B_names, B_Descr, B_TotLen, B_MW, B_SPLen, B_TOTseq] = textread(bBarrelIF,'%s %s %d %d %d %s',-1,'delimiter','\t');
[L_names, L_Descr, L_TotLen, L_MW, L_SPLen, L_TOTseq] = textread(LipoIF,'%s %s %d %d %d %s',-1,'delimiter','\t');

P_MatLen=P_TotLen-P_SPLen;
B_MatLen=B_TotLen-B_SPLen;
L_MatLen=L_TotLen-L_SPLen;

MAX_mat_len=max([P_MatLen;B_MatLen;L_MatLen]);
MAX_SP_len=max([P_SPLen;B_SPLen;L_SPLen]);

%%%%%%%%%%%%%%%   SP LENGTH %%%%%%%%%%%%%%%%

x=12:100:MAX_mat_len
[Pfreq,x_]=hist(P_MatLen,x);
[Bfreq,x_]=hist(B_MatLen,x);
[Lfreq,x_]=hist(L_MatLen,x);

%precent=(freq./Peptides).*100;

Lfreq_perc=(Lfreq./sum(Lfreq)).*100;
Pfreq_perc=(Pfreq./sum(Pfreq)).*100;
Bfreq_perc=(Bfreq./sum(Bfreq)).*100

allfreq=[Bfreq_perc' Pfreq_perc' Lfreq_perc'];

figure(1);
bar3(x',allfreq,'detached');
legend('B-barrel','Periplasmic','Lipoproteins','Location','Best');
zlabel('Percent over all');
ylabel('Mature Length');

%%%%%%%%%%%%%%  MATURE  %%%%%%%%%%%%%%%%%%%%
x=12:4:MAX_SP_len
[Pfreq,x_]=hist(P_SPLen,x);
[Bfreq,x_]=hist(B_SPLen,x);
[Lfreq,x_]=hist(L_SPLen,x);

%precent=(freq./Peptides).*100;

Lfreq_perc=(Lfreq./sum(Lfreq)).*100;
Pfreq_perc=(Pfreq./sum(Pfreq)).*100;
Bfreq_perc=(Bfreq./sum(Bfreq)).*100

allfreq=[Bfreq_perc' Pfreq_perc' Lfreq_perc'];

figure(2);
bar3(x',allfreq,'detached');
legend('B-barrel','Periplasmic','Lipoproteins','Location','Best');
zlabel('Percent over all');
ylabel('SP Length');



end