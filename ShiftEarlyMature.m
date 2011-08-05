function [Combo]=ShiftEarlyMature(SP,MATURE,shift)

mature_len=length(MATURE);
mutations=mature_len-shift;
Combo=cell(mutations,1);

display([SP,MATURE]);
for i=1:1:mutations
    cur_Mutation=[MATURE(shift+1:shift+i),MATURE(1:shift),MATURE(shift+i+1:end)];
    Combo{i}=[SP,cur_Mutation];
%     display(Combo{i});
%     pause;
end

Combo=[{[SP,MATURE]};Combo];

end