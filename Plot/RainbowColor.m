function [RGB_rainbow]=RainbowColor(step,r_g_b)
subDiv=0:step:256;
subDiv_len=length(subDiv);
subDiv_=256:-step:0;
top_=ones(subDiv_len,1)*256;
bottom_=zeros(subDiv_len,1);
Red2Mag=[top_ bottom_ subDiv'];
Mag2Blue=[subDiv_' bottom_ top_];
Blue2Cyan=[bottom_ subDiv' top_];
Cyan2Green=[bottom_ top_ subDiv_'];
Green2Yellow=[subDiv' top_ bottom_];
Yellow2Red=[top_ subDiv_' bottom_];
Colors=[Blue2Cyan(2:end,:);Cyan2Green(2:end,:);Green2Yellow(2:end,:);Yellow2Red(2:end,:);Red2Mag(2:end,:);Mag2Blue(2:end,:)];

if(exist('r_g_b','var')==1)
    switch r_g_b
        case 'b2c'
            Colors=Blue2Cyan(2:end,:);
        case 'c2g'
            Colors=Cyan2Green(2:end,:);
        case 'g2y'
            Colors=Green2Yellow(2:end,:);
        case 'y2r'
            Colors=Yellow2Red(2:end,:);
        case 'r2m'
            Colors=Red2Mag(2:end,:);
        case 'm2b'
            Colors=Mag2Blue(2:end,:);
        case 'y2m'
            Colors=[Yellow2Red(2:end,:);Red2Mag(2:end,:)];
    end
end

RGB_rainbow=Colors./256;

end