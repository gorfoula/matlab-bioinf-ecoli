function [INDEX]=ArraysStrCmp(DatabaseNames,Names)

INDEX=zeros(size(DatabaseNames,1),1);
Proteins=size(Names);

for i=1:1:Proteins
    found=find(strcmp(DatabaseNames,Names{i}));
    if(isempty(found)==0)
        INDEX(found)=1;
    end
end
