function [Catg AAs Catgint AAsint SelSeq Names CatgNames DomProp DomPropIndx DomDimen SPLen MatLen]=FeatureCoding(IFile,type,len,pos,catg,SPlen_distr,aaprop_)

% IFile:        input file
% Ofile:        outputfile
% type:         bin/int/char
% scORnsc:      0 if the input file is with NON Secreted Proteins
%               1 if the input file is with Secreted Proteins
% len:          length of mature Coding
% pos:          starting position on mature
% SPlen_distr:  0/1 not/show SPlen dist
% fasta:        print fasta sequence in fastaout file
% fastaout:     out file for fasta sequence
spyes=sign(pos);
pos=abs(pos);
if(~(pos==0))
    pos=pos-1;
end
len=len+1+pos;  %%% plus one is for Methionine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  READ Dataset %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFile
[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(IFile,'%s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file
[mn, mx, SPlen_max, MAXSP , min_sp , freq]=TotalSPLenDistr(SPLen,'Figures/SP_Len_Dist.bmp',SPlen_distr,'SP Length');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Dataset Coding %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(SPLen(1)>0)  % if the protein is a Secreted one we exclude Meth therefore the SP length is minus one
    SPLen_spare=SPLen-1;
    if(spyes<0)
        len_v=SPLen;
        SPLen=SPLen-SPLen;
    else
        len_v=ones(size(SPLen,1))*len;
        SPLen=SPLen-1;
    end    
else
    len_v=ones(size(SPLen,1))*len;
    SPLen_spare=SPLen;
end
Names=GnNames(TotLen>=SPLen+len);   %Gene Names
SELECTED=SEQ(TotLen>=SPLen+len);    %Select Petides with Seq long enough
samples=length(SELECTED);
TOTALLEN=TotLen(TotLen>=SPLen+len);
SPLen_spare=SPLen_spare(TotLen>=SPLen+len);
SPLen=SPLen(TotLen>=SPLen+len);     %Select corresponding Signal Peptide Lengths
MatLen=TOTALLEN-SPLen;


%%%%% INIT %%%%%%
AAs3D=zeros(samples,20,len-1-pos);
AAs_int=zeros(samples,len-1-pos);
SelSeq=cell(samples,1);
SelSeq_sp=cell(samples,1);
Smast_SP=14;
for i=1:1:samples
    if(spyes==0)
        [AAs3D(i,:,:) SelSeq(i)]=AARepresentation(type,SPLen(i)-Smast_SP+2,SPLen(i)+len-Smast_SP,{SELECTED{i}});  %% Methionine excluded
    else
        [TEMP_REPR SelSeq(i)]=AARepresentation(type,SPLen(i)+2+pos,len_v(i)+SPLen(i),{SELECTED{i}});  %% Methionine excluded
        AAs3D(i,:,1:size(TEMP_REPR,3))=TEMP_REPR;
    end
    
    [AAs_int(i,:) temp]=AARepresentation('int',SPLen(i)+2+pos,len+SPLen(i),{SELECTED{i}});
    [AAs_int_sp SelSeq_sp(i)]=AARepresentation('int',1,SPLen_spare(i),{SELECTED{i}});
end
[CatgNames,Catg3D]=CatgRepr(AAs3D,type,'Features/mature_AAProperties.txt');
[CatgNames_9f,Catg3D_9f]=CatgRepr(AAs3D,type,'Features/AAProperties_new.txt');
[CatgNames,Catg3D_int]=CatgRepr(AAs_int,'int','Features/mature_AAProperties.txt');

[R,C,Z]=size(AAs3D);
[R_catg,C_catg,Z_catg]=size(Catg3D);
[R_catg_9f,C_catg_9f,Z_catg_9f]=size(Catg3D_9f);

AAs3D=reshape(AAs3D,R,C*Z);
Catg3D=reshape(Catg3D,R_catg,C_catg*Z_catg);
Catg3D_9f=reshape(Catg3D_9f,R_catg_9f,C_catg_9f*Z_catg_9f);

Catg=[ones(R,1)*catg Catg3D Catg3D_9f];
AAs=[ones(R,1)*catg AAs3D];
Catgint=[ones(R,1)*catg Catg3D_int];
AAsint=[ones(R,1)*catg AAs_int];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Calculate Net Properties of SP and MAT  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_Charges=8;
max_pho=16;
DomDimen=[];
DomProp=[];
DomPropIndx=[];
if(aaprop_==1)
    if(SPLen(1)==0)
        SPLen=ones(samples,1)*25;
    end
    [Domains AllNames max_Nend max_Hend]=DomainsLoad(Names,'SignalP/SPDomains.txt');
    [Charges Phobic PhoSUm pos_indx pho_indx]=CountPropertiesDomain(SelSeq_sp,AllNames,Names,Domains,max_Nend,max_Hend);
    charge_index=eye(max_Charges+1);
    Pol=charge_index(Charges+1,:);
    pho_index=eye(max_pho+1);
    Pho=pho_index(Phobic+1,:);
    DomProp=[Pol Pho];
    DomPropIndx=[pos_indx pho_indx];
    DomDimen=[max_Charges+1 max_pho+1 max_Nend max_Hend]
end

end