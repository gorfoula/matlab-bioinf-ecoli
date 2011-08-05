function [MW]=MolMass(PeptideSeq,type)
%% INPUT
%%      type:   'm' for monoisotopic
%%              'a'  for average
[OLC TLC FORMULA MW_mono MW_average] = textread('Features/AAMass.txt','%s %s %s %f %f',-1,'delimiter','\t');

aas=size(PeptideSeq,2);
mono=zeros(aas,1);
avg=zeros(aas,1);

for i=1:aas
   found=strcmp(OLC,PeptideSeq(i));
   if(sum(found)>0)
       mono(i)=MW_mono(found);
       avg(i)=MW_average(found);
   end
   
end

switch (type)
    case 'm'
        MW=sum(mono)+18;
    case 'a'
        MW=sum(avg)+18;
    otherwise
        display('Type is: <m> for monoisotopic <a> for average MW');
        return;
end

end