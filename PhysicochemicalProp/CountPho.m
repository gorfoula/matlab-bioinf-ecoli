function [Pho start_idx end_idx matches]=CountPho(SQ,expr,window)
%% INPUT
%%  SQ:     amino acid sequence
%%  expr:   amno acid expression to be matched
%%  window: times of expression repetition found
%% OUTPUT
%%  Pho:    times that expr{window} was found in SQ
%%  start_idx:  table of starting index of matched expr{window}
%%  end_idx:    table of ending index of matched expr{window}
%%  matches:    cell vector of matched strings
if(exist('expr','var')==0)      %if expr is not defined
    expr='[AGLMFIV]';
end

if(exist('window','var')==0)    %if window is not defined then exp matched once
    [start_idx, end_idx, extents, matches]=regexp(SQ,expr);%[FMILVCWATGS]
else
                                                                                                                 %[SGTAWCVLIMF]
    [start_idx, end_idx, extents, matches]=regexp(SQ,[expr,'{',int2str(window(1)),',',int2str(window(end)),'}']);%[FMILVCWATGS]
end
Pho=length(start_idx);

end