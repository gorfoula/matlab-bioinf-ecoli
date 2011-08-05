close all;
S=load('3D secB titration.txt');
timepoints=S(2:end,1);
SecBConc=S(1,2:end);
sampling=0:.5:25;
Units=S(2:end,2:end);
[rows cols]=size(Units);
timepoints_=zeros(length(sampling),cols);
SecBConc_=zeros(length(sampling),cols);
Units_fitted=zeros(length(sampling),SecBConc_);
for i=1:1:cols
%     timepoints_(:,i)=timepoints;
    [x p]=sigmoidfit(timepoints,Units(:,i),sampling);
    plot(sampling,p,'k');
    timepoints_(:,i)=sampling;
    Units_fitted(:,i)=p;

end
for j=1:length(sampling)
    SecBConc_(j,:)=SecBConc;
%     plot(SecBConc,Units(j,:));
%     pause;
end

surfc(timepoints_,SecBConc_,Units_fitted);
xlabel('Time');
ylabel('Concentration');
zlabel('Units');

view(3);