function [ConvTable mask nf]=MaskConv(Matrix,window,sel,angl,norm_)
%% INPUT
%%  Matrix:     vector or array to perform mask convolution
%%              if array mask is applied to each row
%%  window:     length of mask to be applied
%% OUTPUT
%%  ConvTable:  Table resulted from applying the mask

%% INIT
ConvTable=[];
nf=0;

%% Creating Mask
b=randn(10^5,1);
extend=20;
[WeightFactor x]=hist(b,window+extend);
switch sel
    case 'g'
        mask=WeightFactor./10^5;
        mask=mask((extend/2)+1:end-(extend/2));
    case 'l'
        mask=ones(1,window+extend);
        mask=mask((extend/2)+1:end-(extend/2));
    case 'a'
        [mask]=Mask2D(window,angl);
    case 'i'
        [mask]=Mask2D(window,angl);
        mask=double(not(mask));
    otherwise
        display('But selection: <mask type is> is l/g/a');
        return;
end

%% Convolution
conv_=zeros(size(mask,1),length(Matrix)+size(mask,2)-1);
for rot=1:(size(mask,1))
    conv_(rot,:)=conv2(mask(rot,end:-1:1),Matrix);
end
first_conv=((window-1)./2)+1;

if(size(norm_)==1)
    if(norm_(1)==1)
        conv_=abs(conv_)./sum(conv_.^2,2);
        figure;plot(conv_);
    end
    % ConvTable=ConvTable./max(abs(ConvTable));
else
    nf=max(mask*(ones(size(mask,2),1)*norm_(2)));
    conv_=conv_./nf;
end
%%%% area of interest
% ConvTable=conv_(:,window:end);
ConvTable=conv_(:,first_conv:end-first_conv+1);
% ConvTable=conv_(:,1:length(Matrix));
end

function [mask]=Mask2D(window,angl)

aas=1:window;
mask=zeros(180/20,window);
degrees=(aas-1)*100;

for w=1:180/20
    circlepos=mod(degrees,360);
    mask(w,:)=double(circlepos>=angl);
%     mask(w,:)=double(circlepos>=angl)./sum(double(circlepos>=angl));
    degrees=degrees+20;
end



end