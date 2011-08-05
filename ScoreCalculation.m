function [allscores Indx W]=ScoreCalculation(SEQ,SPLen,w,FeatureID,VALUES,len,sel)

MUTATIONS=length(SEQ);
SPLen(SPLen==0)=1;

allscores=zeros(MUTATIONS,1);

W=zeros(1,20*len);
W(FeatureID)=w(1:end-1);
for i=1:1:MUTATIONS    % transfrom to 0-1 coding mode  ALL SP Combinations with specific MATURE
%     display(SEQ{i});
    if(isempty(SPLen))
        [AAs3D SelSeq]=AARepresentation('bin',2,len+1,{SEQ{i}});  %% Methionine excluded
        [Indx SelSeq]=AARepresentation('idx',2,len+1,{SEQ{i}});  %% Methionine excluded
    else
%       [AAs3D SelSeq]=AARepresentation('bin',SPLen(i)+1,len+SPLen(i),{SEQ{i}});  %% Methionine excluded
        [Indx SelSeq]=AARepresentation('idx',SPLen(i)+1,len+SPLen(i),{SEQ{i}});
    end
    
%       AAs_row=reshape(AAs3D,1,VALUES*len);        
    if(sel==20119)
        [CatgNames,Catg3D]=CatgRepr(AAs3D,'bin','Features/mature_AAProperties.txt');
        [CatgNames_9f,Catg3D_9f]=CatgRepr(AAs3D,'bin','Features/AAProperties_new.txt');
        Catg3D=reshape(Catg3D,1,len*11);
        Catg3D_9f=reshape(Catg3D_9f,1,len*9);
        AAs_row=[AAs3D Catg3D Catg3D_9f];
    end
    allscores(i)=sum(W(Indx))+w(end);
%         allscores(i)=AAs_row(:,FeatureID)*w(1:end-1)+w(end);aa=1:1:MUTATIONS;
    end
   
end