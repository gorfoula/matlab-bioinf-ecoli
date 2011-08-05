function [freq]=Freqcalc(array,len)

    cols=size(array,2);
    freq=[];
    for i=1:cols
        [freq_]=hist(array(:,i),len);
        freq_=(freq_./sum(freq_))*100;
        freq=[freq;freq_];
    end

end