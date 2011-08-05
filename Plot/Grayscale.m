function [GREY]=Grayscale(plots)

if(plots>1)
    GREY=0:1/(plots):1;
else
    GREY=0.3;
end
    
    GREY=[GREY' GREY' GREY'];

end