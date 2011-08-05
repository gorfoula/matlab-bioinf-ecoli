Table=ReadTable('JoinFiles/TheoreticalPeptides/CellEnvelopeK12_nocontam_proteinlevel.txt');
Table_OuterMembr=ReadTable('JoinFiles/TheoreticalPeptides/OuterMemb_acc.TXT_SEL.txt');

PIGRAVY=CellTable2Double(Table(2:end,5:6));
PIGRAVY_OuterMembr=CellTable2Double(Table_OuterMembr(2:end,5:6));

ALL={PIGRAVY PIGRAVY_OuterMembr};
HeatMap(ALL,50);
xlabel('GRAVY')
ylabel('PI')