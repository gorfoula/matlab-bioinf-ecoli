function [sum window]=SumOverWindows(indx_logic,window)

sum=0;
for cur_wnd=window
    mask=ones(1,cur_wnd);
    cur_sum=conv(mask,vector);
    m=max(cur_sum);
    pos=find(cur_sum==m,1,'first');
    if(sum<m)
        sum=m;
        diff=pos-cur_wnd+1;
        if(diff>0)
            window=[diff pos];
        else
            window=[1 pos];
        end
    end
end

end
