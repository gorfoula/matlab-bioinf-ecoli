function [names,CatgCors,CATG]=LoadAAFeaturesCatg(file)

fid=fopen(file,'r');

[NumCatg] = textscan(fid,'%d',1,'delimiter','\n');
CATG=NumCatg{1,1};
for i=1:1:NumCatg{1,1}
[names(i)] = textscan(fid,'%s',1,'delimiter','\n');
end

[temp] = textscan(fid,'%d',20,'delimiter',' ');

CatgCors=temp{:,:};
fclose all;

end