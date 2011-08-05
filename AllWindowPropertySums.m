function [MaxWindow MaxDiff]=AllWindowPropertySums(PropertyMatrixA,PropertyMatrixB,windows)

MaxWindow=0;
MaxDiff=0;
for i=windows
    mask=ones(1,i);
    allsumsA=conv2(mask,PropertyMatrixA);
    allsumsB=conv2(mask,PropertyMatrixB);
    [h p ci]=ttest2( allsumsA , allsumsB );
    diff=abs(mean(allsumsA)-mean(allsumsB));
    index=[1:1:length(p)];
    index=index(h==1);
    diff=diff(h==1);
    [curMaxDiff at]=max(diff);
    pos=index(at);
    
    curMaxWindow=[pos-i;pos];
    
    
    if((curMaxDiff>MaxDiff))
        MaxWindow=curMaxWindow;
        MaxDiff=curMaxDiff;        
    end

end