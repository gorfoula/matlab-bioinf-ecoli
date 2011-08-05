function []=FeatureManipulation(SecrIF,CutoIF,type,signal,mature,distr,fasta,fastaout)

%out_26_proteins_SP_EColi_SecA.txt
%out_30_proteins_CP_EColi.txt

AMINOS={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};

%%%%%%%%%%  READ Secreted Sequences  %%%%%%%%%%%%%%%%%%% 
[S_names, S_Descr, S_TotLen, S_MW, S_SPLen, S_TOTseq] = textread(SecrIF,'%s %s %d %d %d %s',-1,'delimiter','\t');
S_MatLen=S_TotLen-S_SPLen;

[mnSM, histmaxSM, SM_at_max, max_SM , min_SM , freqSM]=TotalSPLenDistr(S_MatLen,'Figures/Mature_Len_Dist_Secreted.bmp',0,'Secreted MATURE Length');

% proteins=length(S_names);
% S_SPseq=cell(proteins,1);
% S_MATUREseq=cell(proteins,1);
% for i=1:1:proteins
%     S_SPseq{i}=S_TOTseq{i}(1:S_SPLen(i));
% end
% for i=1:1:length(S_names)
%     S_MATUREseq{i}=S_TOTseq{i}(S_SPLen(i)+1:end);
% end



%%%%%%%%%%  READ CYTO Sequences  %%%%%%%%%%%%%%%%%%%%%%%
[NS_names, NS_Descr, NS_TotLen, NS_MW, NS_TOTseq] = textread(CutoIF,'%s %s %d %d %s',-1,'delimiter','\t');
[mnNSM, histmaxNSM, NSM_at_max, max_NSM , min_NSM , freqNSM]=TotalSPLenDistr(NS_TotLen,'Figures/Mature_Len_Dist_NSecreted.bmp',0,'Non Secreted MATURE Length');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s=fopen(['All_Secreted_Proteins_Data.txt'],'w');
% fprintf(s,'Name\tDescription\tTotalLen\tSPLen\tMATLen\tMW\tSP_Sequence\tMAT_Sequence\n');
% for i=1:1:proteins
%     fprintf(s,'%s\t%s\t%d\t%d\t%d\t',S_names{i},S_Descr{i},S_TotLen(i),S_SPLen(i),S_TotLen(i)-S_SPLen(i));
%     fprintf(s,'%d\t%s\t%s\n',S_MW(i),S_SPseq{i},S_MATUREseq{i});
% end
% 
% fclose(s);
% fclose all;


%[mn, mx, SPlen_max, max_sp , min_sp , freq]
[mn, max, SPlen_max, max_sp , min_sp , freq]=TotalSPLenDistr(S_SPLen,'Figures/SP_Len_Dist.bmp',0,'SP Length');

if(signal==0)
    SP_len=SPlen_max;
elseif (freq(signal)==0)
    display(['Peptides with SP length  <',int2str(signal),'> does not exist.']);
    display(['MAX: >',int2str(max_sp),'  MIN: >',int2str(min_sp)]);
    return;
else
    SP_len=signal;
end
MAT_len=mature;

[MAT_len SP_len]

len=MAT_len+max_sp+1;%min(TotLen-SP_Len);%min(SP_Len);
pos=3;
len=103+pos;

NON_SECRETED=NS_TOTseq(NS_TotLen>=len);
NON_SECRETED_Names=NS_names(NS_TotLen>=len);
SelSeqCyto=cell(length(NON_SECRETED),1);
[AAs_Cytopl SelSeqCyto]=AARepresentation(type,2+pos,len,NON_SECRETED);

% SECRETED_Names=S_names(S_TotLen>=len);
% SECRETED=S_TOTseq(S_TotLen>=len);
% [AAs_Secr SelSeqPeri]=AARepresentation(type,2,len,SECRETED);
% [AAs_Secr_int SelSeqPeri_int]=AARepresentation('int',2,len,SECRETED);

SECRETED_Names=S_names(S_TotLen>=S_SPLen+len-2);
SECRETED=S_TOTseq(S_TotLen>=S_SPLen+len-2);
samples=length(SECRETED);
S_SPLen=S_SPLen(S_TotLen>=S_SPLen+len-2);
AAs_Secr=zeros(samples,20,len-2-pos);
AAs_Secr_int=zeros(samples,len-2-pos);
SelSeqPeri=cell(samples,1);
for i=1:1:samples
    [AAs_Secr(i,:,:) SelSeqPeri(i)]=AARepresentation(type,S_SPLen(i)+1+pos,len+S_SPLen(i)-1,{SECRETED{i}});
    [AAs_Secr_int(i,:) SelSeqPeri(i)]=AARepresentation('int',S_SPLen(i)+1+pos,len+S_SPLen(i)-1,{SECRETED{i}});    
end
[R,C,Z]=size(AAs_Secr);

[Names_S,Catg_Secr]=CatgRepr(AAs_Secr,type,'Features/AAProperties10.txt');
[Names_S_int,Catg_Secr_int]=CatgRepr(AAs_Secr_int,'int','Features/AAProperties10.txt');
[Names_C,Catg_Cytopl]=CatgRepr(AAs_Cytopl,type,'Features/AAProperties10.txt');

[R,C,Z]=size(AAs_Secr);
[R_catg,C_catg,Z_catg]=size(Catg_Secr);
CATG=length(Names_C);

if (strcmp(type,'bin'))
    CYTOPLASMIC=length(AAs_Cytopl(:,1,1));
    AAs_Secr=reshape(AAs_Secr,R,C*Z);
    AAs_Cytopl=reshape(AAs_Cytopl,CYTOPLASMIC,C*Z);
    Catg_Secr=reshape(Catg_Secr,R_catg,C_catg*Z_catg);
    Catg_Cytopl=reshape(Catg_Cytopl,CYTOPLASMIC,C_catg*Z_catg);
    cat_c_indx=1:len*C_catg;
    aas_c_indx=1:len*C;
elseif(strcmp(type,'int'))
    cat_c_indx=(1:len);
    aas_c_indx=(1:len);
    
end

[R,C,Z]=size(AAs_Secr);

Catg_Secr_sp=[ones(R,1) Catg_Secr];
Catg_Cytopl_sp=[zeros(CYTOPLASMIC,1) Catg_Cytopl];
AAs_Secr_sp=[ones(R,1) AAs_Secr];
AAs_Cytopl_sp=[zeros(CYTOPLASMIC,1) AAs_Cytopl];

%%%%%%%%  WRITE in file Sequence Representation  %%%%%%%%
Tot_6f_sp=[Catg_Secr_sp;Catg_Cytopl_sp];
Tot_20f_sp=[AAs_Secr_sp;AAs_Cytopl_sp];

Initial=[CYTOPLASMIC+R ((20*Z)+1) 2]
Initial_cat=[CYTOPLASMIC+R (CATG*Z)+1 2]

dlmwrite(['Gems/_P',SecrIF,'_t',type,'_.txt'], Initial,'delimiter','\t');  %  write the Secreted proteins
dlmwrite(['Gems/_P',SecrIF,'_t',type,'_.txt'], Tot_20f_sp,'delimiter','\t','-append');  %  write the Secreted proteins

dlmwrite(['Gems/',int2str(CATG),'f_P',SecrIF,'_t',type,'_.txt'], Initial_cat,'delimiter','\t');  %  write the Secreted proteins
dlmwrite(['Gems/',int2str(CATG),'f_P',SecrIF,'_t',type,'_.txt'], Tot_6f_sp,'delimiter','\t','-append');  %  write the Secreted proteins

%%%%%%%%%%%%%%  FASTA LOGO  %%%%%%%%%%%%%%%%%%%%%
if(fasta==1)
    Selected=[SelSeqCyto];
    FastaOut(['fASTA/',fastaout],[NON_SECRETED_Names],Selected,1,len-2-pos);
%     Selected=[SelSeqPeri;SelSeqCyto];
%     FastaOut(['fASTA/',fastaout],[SECRETED_Names;NON_SECRETED_Names],Selected,1,len-2-pos);
end
%%%%%%%%%%%  Distribution  %%%%%%%%%%%%%%%%%


if(distr==1)    
    AminoDistrFig(Catg_Secr_int,Names_S,1,['Figures/',int2str(length(Names_S)),'Catg_SP',int2str(SP_len),'_']);
    AminoDistrFig(AAs_Secr_int,AMINOS,1,['Figures/Catg_SP',int2str(SP_len),'_']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all;

end