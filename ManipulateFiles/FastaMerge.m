function []=FastaMerge(IFasta,IFasta2)
% FastaSubSelection('Data/InnerMembrane_SPID.txt','Data/Malvina/InnerMembraneEcoli.fasta','Data/Proteins_IDs.txt')

found=regexp(IFasta,'[/]');
dirfasta=IFasta(1:found(end));

[HEAD SEQUENCE]=FastaRead(IFasta);
[HEAD2 SEQUENCE2]=FastaRead(IFasta2);

FastaWrite([dirfasta,'Merged.fasta'],HEAD,SEQUENCE,70,'w');
PROTEINS2=length(HEAD2);
HEAD_NEW=cell(1,PROTEINS2);
SEQUENCE_NEW=cell(1,PROTEINS2);

[MNEM_st MNEM_en MNEM_ext MNEM]=regexp(HEAD,'[a-zA-Z0-9]{3,5}_(ECOLI)'); % Gene names e.g. GN=phoA
for i=1:1:length(MNEM)
    if(isempty(MNEM{i})==0)
        MNEM{i}=MNEM{i}{1};
    end
end
[MNEM2_st MNEM2_en MNEM2_ext MNEM2]=regexp(HEAD2,'[a-zA-Z0-9]{3,5}_(ECOLI)'); % Gene names e.g. GN=phoA
for i=1:1:length(MNEM2)
    if(isempty(MNEM2{i})==0)
        MNEM2{i}=MNEM2{i}{1};
    end
end


for i=1:1:PROTEINS2
    if(isempty(MNEM2{i}))
        HEAD_NEW{i}=HEAD2{i};
        SEQUENCE_NEW{i}=SEQUENCE2{i};
    else
        indx_seq=find(strcmpi(MNEM,MNEM2{i}));
    end
    if(isempty(indx_seq))
        HEAD_NEW{i}=HEAD2{i};
        SEQUENCE_NEW{i}=SEQUENCE2{i};
    end
end

FastaWrite([dirfasta,'Merged.fasta'],HEAD_NEW,SEQUENCE_NEW,70,'a');

end