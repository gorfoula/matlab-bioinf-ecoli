function [positions sequences]=IsTherePhoPatch(Name,SEQ,windows)
%% INPUT
%%  Name:       Name cell vevtor of sequences
%%  SEQ:        sequence vector
%%  windows:    length window of hydrophobic patches
%% OUTPUT
%%  positions:   double array => length of patch | start of patch | end of patch
%%  sequences:   cell array => Name | matched sequences
[Pho start_indx end_indx sequences]=CountPho(SEQ,'[LIVFCMA]',windows);
positions=[(end_indx'-start_indx'+1) start_indx' end_indx'];
sequences={Name,sequences};

end