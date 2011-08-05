function [FEATURES_color]=AssignFeatureColor(FEATURES)

FEATURES_color=zeros(FEATURES,3);

if(FEATURES>20 && FEATURES<=26)
    for i=21:1:FEATURES
        FEATURES_color(i,:)=([200 128 0]./256);   % brown
    end
elseif(FEATURES>26)
    for i=41:1:FEATURES
        FEATURES_color(i,:)=([200 128 0]./256);   % brown
    end
    %{'Bulk' 'PosCh' 'NegCh' 'Pho' 'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S'
    %'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E' 'h' 's' 'R' 'K' 'D' 'E' 'x' 'q' 'H' 'M' 'W' '@' '+' 's' 'o' 'h' 'b' 'P' 'q' 'C'};
    FEATURES_color(21,:)=[0 1 0];    %V,L,I,M,F green
    FEATURES_color(22,:)=[1 1 0];
    FEATURES_color(23,:)=[0 0 1];
    FEATURES_color(24,:)=[0 0 1];
    FEATURES_color(25,:)=[1 0 0];
    FEATURES_color(26,:)=[1 0 0];
    FEATURES_color(27,:)=[0 0 0];
    FEATURES_color(28,:)=([150 144 144]./255);    %q  black
    FEATURES_color(29,:)=([255 128 0]./256);
    FEATURES_color(30,:)=([200 128 0]./256);
    FEATURES_color(31,:)=[1 0 1];
    
    FEATURES_color(32,:)=[1 0 0];   %D,E red
    FEATURES_color(33,:)=[0 0 1];   %K,R blue
    FEATURES_color(34,:)=[1 1 0];   %small yellow
    FEATURES_color(36,:)=[0.7 1 0];    %h green
    FEATURES_color(37,:)=[1 0 1];    % bulky magenta
    FEATURES_color(38,:)=([255 128 0]./256);    % orange P
    FEATURES_color(39,:)=([150 144 144]./255);    %N,Q,H  grey 
end


FEATURES_color(19,:)=[1 0 0];   %D,E red
FEATURES_color(20,:)=[1 0 0];

FEATURES_color(16,:)=[0 0 1];   %K,R blue
FEATURES_color(17,:)=[0 0 1];

FEATURES_color(2,:)=[1 1 0];   %A,G yellow
FEATURES_color(9,:)=[1 1 0];

% FEATURES_color(10,:)=[0 0 0];    %T,S  black
% FEATURES_color(15,:)=[0 0 0];
% FEATURES_color(15,:)=[0 0 0];

FEATURES_color(1,:)=[0 1 0];    %V,L,I,M,F green
FEATURES_color(3,:)=[0 1 0];
FEATURES_color(4,:)=[0 1 0];
FEATURES_color(5,:)=[0 1 0];
FEATURES_color(7,:)=[0 1 0];

FEATURES_color(8,:)=[1 0 1];    % W,Y magenta
FEATURES_color(14,:)=[1 0 1];

FEATURES_color(6,:)=([255 128 0]./256);    %P  orange

FEATURES_color(12,:)=([150 144 144]./255);    %N,Q,H  grey
FEATURES_color(13,:)=([150 144 144]./255);    
FEATURES_color(18,:)=([150 144 144]./255);


end