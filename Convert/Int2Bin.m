function [binary]=Int2Bin(Integer)

binary=[];
multiplier=Integer;
    while(sum(multiplier)>0)
        residuum=rem(multiplier,2);
        multiplier=floor(multiplier./2);
        binary=[logical(residuum) logical(binary)];
        
    end
end