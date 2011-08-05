function []=AllCategoriesCoding(IFile,Ofile,type,len,pos,SPlen_distr,distr,fasta,fastaout,DomProp_flag)

%%%%% READ FILE NAMES  %%%%%%%%%%
[FileNames] = textread(IFile,'%s',-1,'delimiter','\t');

founddot=regexp(IFile,'[.]');
filename=IFile(founddot(end)+1:end)

AMINOS={'M' 'A' 'V' 'L' 'I' 'P' 'F' 'W' 'G' 'S' 'C' 'N' 'Q' 'Y' 'T' 'K' 'R' 'H' 'D' 'E'};

for  i=1:1:length(FileNames)
   [Catg_t AAs_t Catg_int_t AAs_int_t SelSeq_t Names_t Catg_Names DomProp_t DomPropIndx_t DomDimen SPLen_t MatLen_t]=FeatureCoding(FileNames{i},type,len,pos,i-1,SPlen_distr,DomProp_flag);
   samples=length(Names_t);
   {FileNames{i}}
   if(i==1)
       Catg_secr=Catg_t;
       AAs_secr=AAs_t;
       Catg_secr_int=Catg_int_t;
       AAs_secr_int=AAs_int_t;

       SelSeq_secr=SelSeq_t;
       Names_secr=Names_t;
      
       DomProp=[ones(samples,1)*i DomProp_t];
       DomPropIndx=[ones(samples,1)*i DomPropIndx_t];
       
       SPLen=SPLen_t;
       MatLen=MatLen_t;
       
   else
       Catg_secr=[Catg_secr;Catg_t];
       AAs_secr=[AAs_secr;AAs_t];
       Catg_secr_int=[Catg_secr_int;Catg_int_t];
       AAs_secr_int=[AAs_secr_int;AAs_int_t];
       
       SelSeq_secr=[SelSeq_secr;SelSeq_t];
       Names_secr=[Names_secr;Names_t];
       
       DomProp=[DomProp;[ones(samples,1)*i DomProp_t]];
       DomPropIndx=[DomPropIndx;[ones(samples,1)*i DomPropIndx_t]];
       
       SPLen=[SPLen;SPLen_t];
       MatLen=[MatLen;MatLen_t];
   end
 
end

% PropertyDist(1,BulkSum,'Bulkiness','SP Bulkiness compared to Cytoplasmic Proteins (32 first aa)')
if(DomProp_flag==1)
    Tot_20=[AAs_secr(:,1) DomProp(:,2:end) DomPropIndx(:,2:end)];
elseif(DomProp_flag==20119)
    Tot_20=[AAs_secr(:,1) AAs_secr(:,2:end) Catg_secr(:,2:end)];
else
    Tot_20=[AAs_secr(:,1) AAs_secr(:,2:end)];
end
% for i=1:1:6
% iter=10^7;
% [FS]=FeatureSelection(Tot_20,[0 1],[0 1],iter,['Data/',Ofile,'_jointdist.txt'])
% found=regexp(Ofile,'[.]');
% index=1:1:FS.totalfeatures;
% dlmwrite(['Gems Results/NewSplit/',Ofile(1:found-1),'_.txt'],[index' (find(FS.state==1)+1)'],'delimiter','\t');
% SelectedFeatures(['Gems Results/NewSplit/',Ofile(1:found-1),'_.txt'],0,'Gems/CytoVSsecreted.gemsin','Gems/CytoVSsecreted.test',1);
% end


[R,C_aas]=size(Tot_20);
% [R,C_catg]=size(Tot_f);

Initial=[length(Tot_20(:,1)) C_aas length(FileNames)];
% Initial_cat=[length(Tot_f(:,1)) C_catg length(FileNames)];
found=regexp(Ofile,'[.]');
if(DomProp_flag==1)
    dlmwrite([Ofile(1:found-1),'(',num2str(DomDimen(1)),'.',num2str(DomDimen(2)),'.',num2str(DomDimen(3)),'.',num2str(DomDimen(4)),')',Ofile(found:end)], Initial,'delimiter','\t');  %  write the Secreted proteins  NumPos/NumPho/indxPos/indxPho
    dlmwrite([Ofile(1:found-1),'(',num2str(DomDimen(1)),'.',num2str(DomDimen(2)),'.',num2str(DomDimen(3)),'.',num2str(DomDimen(4)),')',Ofile(found:end)], Tot_20,'delimiter','\t','-append');  %  write the Secreted proteins
else
    dlmwrite(Ofile, Initial,'delimiter','\t');  %  write the Secreted proteins
    dlmwrite(Ofile, Tot_20,'delimiter','\t','-append');  %  write the Secreted proteins
%     dlmwrite(['Gems/',filename,'.test'], Initial_cat,'delimiter','\t');  %  write the Secreted proteins
%     dlmwrite(['Gems/',filename,'.test'], Tot_f,'delimiter','\t','-append');  % write the Secreted proteins
end

%%%%%%%%%%%%%%  FASTA LOGO  %%%%%%%%%%%%%%%%%%%%%
if(fasta==1)
    Selected=[SelSeq_secr];
    FastaLine(fastaout,[Names_secr],Selected,1,len);
end
%%%%%%%%%%%  Distribution  %%%%%%%%%%%%%%%%%


if(distr==1)
    SECR=Catg_secr_int(Catg_secr_int(:,1)==1,2:end);
    CYTO=Catg_secr_int(Catg_secr_int(:,1)==0,2:end);
    AminoDistrFig(CYTO,Catg_Names,1,['Figures/',int2str(length(Catg_Names)),'Catg_SP',int2str(len),'_']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose all;

end