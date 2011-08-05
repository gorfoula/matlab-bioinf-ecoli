clear all;
filename='SecABact_3200';
[HEAD SEQUENCE]=FastaRead(['Data/SecA/',filename,'.fasta']);
for i=1:1:500
    if(length(SEQUENCE{i})>=16)
        FastaWrite(['Data/SecA/',filename,'_500_15.fasta'],{HEAD{i}},{SEQUENCE{i}(1:16)},70,'a');
    end
end
fclose all;