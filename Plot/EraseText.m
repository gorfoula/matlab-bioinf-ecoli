function []=EraseText(FH)

for j=1:length(FH)
    delete(FH(j));
end

end