function [Combo]=ChopEarlyMature(SP,MATURE,step,len)

sp_len=length(SP);
mature_len=length(MATURE);
mutations=floor((mature_len-len)./step);
Combo=cell(mutations,1);

display([SP,MATURE]);
for i=1:1:mutations
    cur_Mature=MATURE(i*step+1:i*step+len-sp_len);
    Combo{i}=[SP,cur_Mature];
%     display(Combo{i});
%     pause;
end

Combo=[{[SP,MATURE(1:len-sp_len)]};Combo];
end