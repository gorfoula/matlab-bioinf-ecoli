function [m2z PeptMass PossibleCharges]=Mass2Charge(Sequence)

N=size(Sequence,1);
m2z=zeros(N,1);
PeptMass=zeros(N,1);
PossibleCharges=zeros(N,1);

for i=1:1:N
    [matched]=regexp(Sequence{i},'[HKR]');
    PossibleCharges(i)=size(matched,2)+1;
    PeptMass(i)=MolMass(Sequence{i},'m');
    m2z(i)=PeptMass(i)/PossibleCharges(i);
end


end