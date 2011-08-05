function [disorderRegions count totalDisRegns]=ReadDisorderAreas(DataFile,DisorderFile)

[GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(DataFile,'%s %s %d %f %d %s',-1,'delimiter','\t');  %%% read dataset file
[lines] = textread(DisorderFile,'%s',-1,'delimiter','\t');
ind=0;count=0;spend=0;
disorderRegions=[];
totalDisRegns=[];
peptide_cur=0;
for i=1:1:length(lines)
    found=regexp(lines{i},'[>\n]');
    region=regexp(lines{i},'[1-9]+[-][1-9]+','once');
    if(isempty(found)==0)
       ind=find(strcmpi(GnNames,lines{i}(found(1)+1:end)));
       GnNames{ind};
       spend=SPLen(ind);
       if(peptide_cur>0)
            totalDisRegns(peptide_cur)=cur_tot_regns;
       end
       peptide_cur=peptide_cur+1;
       cur_tot_regns=0;
    end
    if(isempty(region)==0)
        separator=regexp(lines{i},'-');
        start_=str2num(lines{i}(1:separator-1))-spend;
        end_=str2num(lines{i}(separator+1:end-1))-spend;
        if(start_>=1 && start_<=200)
            count=count+1;
            disorderRegions(count,1:3)=[end_-start_+1 start_ end_];
            cur_tot_regns=cur_tot_regns+1;
        end
    end
end
end