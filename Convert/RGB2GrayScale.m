function [GREY]=RGB2GrayScale(RGB)
R=RGB(:,1);
G=RGB(:,2);
B=RGB(:,3);

GREY=0.3*R + 0.59*G + 0.11*B;

GREY=[GREY GREY GREY];


end