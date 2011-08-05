function [SLSQ]=SequenceSplit(SQ,splits)

len=length(SQ);
lines=ceil(len/splits);
SLSQ=cell(1,lines);

for i=1:lines
    if(i==lines)
        SLSQ{i}=SQ((i-1)*splits+1:end);
    else
        SLSQ{i}=SQ((i-1)*splits+1:(i-1)*splits+splits);
    end
end

end