function [mn, mx, SPlen_max, max_sp , min_sp , freq]=TotalSPLenDistr(SP_Len,fileFig,show,Xlabel)

[Peptides Col]=size(SP_Len);
mn=mean(SP_Len);
max_sp=max(SP_Len);
min_sp=min(SP_Len);
%%%%%%%%%% hist SP_length  %%%%%%%%%%%%%%%%%
step=round((max_sp-min_sp)*.05);
step=1;
x=min_sp:step:max_sp;
[freq,x]=hist(SP_Len,x);
mx=max(freq);  % max value
SPlen_max=x(freq==max(freq))  %SP length with more peptides

if(show==1)
f_handle=figure(1);
precent=(freq./Peptides).*100;
bar(x,precent);  %%%  precentage of all proteins having this property/amino acid
title(['MEAN: ',num2str(mn),' MAX: ',int2str(SPlen_max),' std: ',int2str(std(SP_Len)),' sqrt(var): ',int2str( sqrt( var(SP_Len) ) )]);
xlabel(Xlabel);
ylabel('Precentage over all');
hold on;
x_=[0:.1:max_sp];
plot(mean(SP_Len)*ones(1,length(x_)),x_,'g.');
ylim([0 max(precent)]);
xlim([min_sp-1 max_sp+1]);
hold off;
saveas(f_handle,fileFig,'bmp');
pause;
end


end