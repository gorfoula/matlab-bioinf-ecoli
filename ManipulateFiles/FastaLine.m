function []=FastaLine(file,Names,SP,startpos,endpos)

file_handle=fopen(file,'w');
[R,C,Z]=size(SP);
for i=1:1:R
    if(length(SP{i,1})>=endpos)
        fprintf(file_handle,'>%s\n%s\n', Names{i,1},SP{i,1}(startpos:endpos));
    else
        fprintf(file_handle,'>%s\n%s\n', Names{i,1},SP{i,1}(startpos:end));
    end
end

fclose(file_handle);

end