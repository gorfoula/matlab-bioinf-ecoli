function [allsums]=WindowPropertySum(PropertyMatrix,window)

    mask=ones(1,window);
    allsums=conv2(mask,PropertyMatrix);
    
end
