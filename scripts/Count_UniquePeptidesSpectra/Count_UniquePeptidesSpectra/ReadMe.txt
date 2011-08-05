%% This function given the output of scafold it calculates
%%      -> Unique Peptides
%%      -> Unique Spectra
%%      -> Total Number of spectra
%%
%% INPUT
%%      file:   path of output file of scafold, in tab delimited version
%%      filedb: database file
%%      index:  [H X Y Z1 Z2 A B C]
%%              H:  how many lines are occupied by header
%%              X:  column number of Protein Accession (1)
%%              Y:  column number of Peptide sequence (2)
%%              Z1: Protein MW (3)
%%              Z2: column number of Peptide MW (4)
%%              A:  column number 1H counts (5)
%%              B:  column number 2H counts (6)
%%              C:  column number 3H counts (7)
%%
%%  EXAMPLE: CountSpectra('Peptide Report for PRO008_1DG010_L1_BL21_allmods.txt','TableS4_BL21_DE3_ProteomeV8.txt');