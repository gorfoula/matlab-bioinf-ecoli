function [data]=ReadSelectedFeatures(DataINfile,FeatureID,NUMofFEATURES)

fid=fopen(DataINfile,'r');
[PARAM] = textscan(fid,'%d',3,'delimiter','\t');
PEPTIDES=PARAM{1}(1);
FEATURES=PARAM{1}(2);

data=zeros(PEPTIDES,FEATURES);
for i=1:1:PEPTIDES
temp= textscan(fid,'%d',FEATURES,'delimiter','\t');
data(i,:)=temp{1};
end


end