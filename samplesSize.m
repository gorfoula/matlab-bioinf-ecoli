function []=samplesSize(E)

z=1.96;  %%% 95%sure
s=0.89; %%%% standard deviation of score

if(exist('E','var')==0)
    E=0.1:0.1:0.5;
end
    
n=((z*s)./E').^2;

table=[n E'];
dlmwrite(['Data/Sample size.txt'], table,'delimiter','\t');  %  write the Secreted proteins

end