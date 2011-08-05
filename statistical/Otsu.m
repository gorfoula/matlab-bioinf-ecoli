function [Threshold,index]=Otsu(h,pixels,par)

bins=50;

[histogram,x]=hist(h,bins);        

tmp=max(max(x));
h_=h(h<tmp);
Threshold=max(max(h_));
index=(h>Threshold);

prev_between_variance=0;
minh=min(min(h));


while(Threshold>=(minh+.0001)) 
    index=(h>Threshold);% & h<200;     %% Set threshold
    background=not(index);

%     mean_fg=mean(h(index));
%     mean_bg=mean(h(background));

    Num_fg_pixels=sum(histogram(x>Threshold)/pixels);
    Num_bg_pixels=1-Num_fg_pixels;
    
    mean_fg=(x(x>Threshold)*histogram(x>Threshold)')/(pixels*Num_fg_pixels);
    mean_bg=(x(x<=Threshold)*histogram(x<=Threshold)')/(pixels*Num_bg_pixels);
    
    if(pixels*Num_fg_pixels==0)
        pixels
        Num_fg_pixels
        histogram(x>Threshold)
        histogram(x>Threshold)/pixels
    end
        


    % between_var=no*nb(mb-mo)^2
    % no: number of foreground pixels
    % nb: number of background pixels
    % mb,mo: mean values of bg and fg pixels
    dif=mean_bg-mean_fg;
    between_variance=Num_fg_pixels*Num_bg_pixels*(dif^2);

    prev_between_variance;     
        
    if (between_variance<prev_between_variance)
            break;
    end 
        prev_between_variance=between_variance;
        if(par=='h')
            Threshold=Threshold-.01;
        else
            Threshold=max(x(x<Threshold));
        end
    
    
end