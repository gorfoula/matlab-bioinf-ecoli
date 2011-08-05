function []=SPLengthHistogram(IFile)

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file

[mn, mx, SPlen_max, MAXSP , min_sp , freq]=TotalSPLenDistr(SPLen,'Figures/SP_Len_Dist.bmp',1,'SP Length');


end