function [indx]=Classify_Kmeans(S,ctrs)

if(size(S,2)~=size(ctrs,2))
    display('Number of parameters of X and ctrs must be the same');
    return;
end

samples=size(S,1);
classes=size(ctrs,1);
params=size(ctrs,2);
distance=zeros(samples,classes);
for i=1:1:classes
    [C]=meshgrid(ctrs(i,:),1:1:samples);
    distance(:,i)=sqrt(sum((S-C).^2,2));
end

min_dist=(distance==(meshgrid(min(distance'),1:1:classes)'));
[class_indx]=meshgrid(1:1:classes,1:1:samples);

indx=sum(class_indx.*min_dist,2);
end