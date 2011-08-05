function [Domains genename Catg]=ExtractDomainsSignalP(lines,catg)

NumLines=length(lines);

res_flag=0;
count=1;
table=zeros(70,3);
peptides=0;
tableend=[];
for i=1:1:NumLines
    
    name=regexp(lines{i},'[>]');
    results=regexp(lines{i},'[#]','once');
    probcol=regexp(lines{i},'[0-9]+.[0-9]+');
    
%     lines{i}
    if(not(isempty(name)))
      if(res_flag && tableflag)
          GN
        if(strcmp(GN,'gspI'))
            display('Hi!');
        end
        
        prob_thr=table(:,3)>.3;
        c_cur_start=find(prob_thr==1,1,'first');
        c_cur=find(prob_thr==1,1,'last');
        c_cur_true=find(prob_thr(c_cur_start:c_cur)==0,1,'first');
        if(isempty(c_cur_true)==0)
            c_cur=c_cur_true+c_cur_start-2;
        end
        table=table(1:c_cur,:);
        
        [m index]=max(table');
        n_end=find(index==2,1,'first')-1;
        h_end=find(index==2,1,'last');
        c_end=find(index==3,1,'last');
        Domains(peptides,1:6)=[n_end h_end c_end n_end (h_end-n_end) (c_end-h_end)];
        count=1;res_flag=0;tableflag=0;
      end
      GN=lines{i}(name+1:end);
    end
    if(not(isempty(results)))
      res_flag=1;
      peptides=peptides+1;
      genename{peptides}=GN;
      Catg{peptides}=catg;      
    end
    if(res_flag && not(isempty(probcol)))
%         lines{i}
      ndomain=str2double(lines{i}(probcol(3):probcol(4)-1));
      hdomain=str2double(lines{i}(probcol(4):probcol(5)-1));
      cdomain=str2double(lines{i}(probcol(5):end));
      table(count,:)=[ndomain hdomain cdomain];
      tableflag=1;
      count=count+1;
%       tableend=regexp(lines{i}(1:end-1),['>',genename{peptides}],'once');
    end
end
genename=genename';
Catg=Catg';
end