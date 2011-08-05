function [filename directory]=IsolateFileName(IFPath)

NoNames=size(IFPath,1);
filename=cell(NoNames,1);
directory=cell(NoNames,1);

for f=1:1:NoNames
    
[dir]=regexp(IFPath{f},'[/]');               %%% directory of files
[dot]=regexp(IFPath{f},'[.]');              %%% suffix of file
filename{f}=IFPath{f};
if (isempty(dir)==0 && isempty(dot)==0)
    filename{f}=IFPath{f}(dir(end)+1:dot(end)-1);
    directory{f}=IFPath{f}(1:dir(end));
end

end
end