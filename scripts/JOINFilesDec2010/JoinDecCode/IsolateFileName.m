function [filename directory]=IsolateFileName(IFPath)

[dir]=regexp(IFPath,'[/]');               %%% directory of files
[dot]=regexp(IFPath,'[.]');              %%% suffix of file
filename=IFPath;
if (isempty(dir)==0 && isempty(dot)==0)
    filename=IFPath(dir(end)+1:dot(end)-1);
    directory=IFPath(1:dir(end));
end

end