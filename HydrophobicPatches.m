function [positions sequences patches AllPosPho]=HydrophobicPatches(DatasetF,len,windows,type)
%% INPUT
%%  DatasetIF:  Dataset in file, file with specific conformation
%%  windows:    lenght range of hydropho patches
%% OUTPUT
%%  positions:  array => length of patch | start of patch | end of patch
%%  sequences:  cell array => Name | matched sequences
%%  patches:    Number of patches found[text_] = textread(IFNames{1},'%s',-1,'delimiter','\n'); % read all lines
%% Read Database File
[text_] = textread(DatasetF,'%s',-1,'delimiter','\n'); % read all lines
header_long=text_{1};
[start_idx, end_idx, extents, matches, tokens, names, Table_] = regexp(text_,'[\t]'); % split coloumns
Table_Secr=CellTable2StrTable(Table_); % convert Cell table to string table

GnNames=Table_Secr(:,1);
TotLen=CellTable2Double(Table_Secr(:,3));
SPLen=CellTable2Double(Table_Secr(:,5));
SEQ=Table_Secr(:,6);

% [GnNames, Descr, TotLen, MW, SPLen, SEQ] = textread(DatasetF,'%s %s %d %s %d %s',-1,'delimiter','\t');  %%% read dataset file

if(SPLen(1)==0)
    SPLen=SPLen+1;
end
SEQ=SEQ(TotLen>(len+SPLen));
Names=GnNames(TotLen>(len+SPLen));
SPLen=SPLen(TotLen>(len+SPLen));
Peptides=length(SEQ);
positions=[];AllPosPho=[];Profile=[];
AllPosPho_cur=[];sequences_cur=[];
patches=zeros(Peptides,1);
sequences={};

for i=1:1:Peptides
    mature=[SPLen(i)+1 SPLen(i)+len];
    cur_seq=SEQ{i}(mature(1):mature(2));
    switch (type)
        case 'seq'
            %% Sequences of continous hydrophobic amino acids
            [positions_cur sequences_cur]=IsTherePhoPatch(Names{i},cur_seq,windows);
        case 'scale'
            ConvTable=zeros(1,len);
            %% Pho Scale
            [Profile_cur]=ScaleProfile('K',Names{i},cur_seq,windows(3),0);

            %% find areas
%             Profile_cur=(Profile_cur)-.8;
%             PlotProfile(1,Names{i},cur_seq,Profile_cur);pause;
            
            pos_area=find(diff([0 0 double(Profile_cur>0) 0 0],2)>0);  %% 2h paragwgos = shmeio kampis (allagh proshmou allazei kampylothta)
%             figure(1); hold on; xlim([-1 102]);stem(0:103,pos_area(1:104));hold off;
            W=pos_area(2:end)-pos_area(1:end-1)-1;
            width=W(1:2:length(W));
            start_pos=pos_area(1:2:length(W))+1;
            center=start_pos+floor(width./2);
            
            width_limits=(width>=windows(1) & width<=windows(2));
            center=center(width_limits);
            patch_start=start_pos(width_limits);
            width=width(width_limits);
            positions_cur=[width' patch_start' (patch_start+width-1)' center'];         % table => center of patch|length of patch           

        otherwise
            display('=>Wrong type Hydrophobic patches scale/seq <HydrophobicPatches.m>');
            return;
    end

    if(isempty(positions_cur)==0)
        [AllPosPho_cur]=AllPositionsArray(positions_cur); % double vectror with all positions whith pho Amino acid
    end
    %% append current sequence calculations
    if(isempty(positions_cur)==0)
        AllPosPho=[AllPosPho AllPosPho_cur];
        positions=[positions;positions_cur];
        sequences=[sequences;sequences_cur];
        patches(i)=length(positions_cur(:,1));
    end
    
end

end

function [width center]=FindAreas(limits,Profile)

areas=length(limits);
if(mod(areas,2))
    limits=[limits length(Profile)];
end
areas=(areas+mod(areas,2))./2;

width=zeros(1,areas);
center=zeros(1,areas);

pos=1;
for i=1:(areas*2)-1
    width(pos)=max(Profile(limits(i):limits(i+1)));
    center(pos)=limits(i)+floor((limits(i+1)-limits(i))/2);
    pos=pos+1;
end

end

function [AllPosPho]=AllPositionsArray(positions)
%% INPUT
%%  positions:   double array => length of patch | start of patch | end of
%%  patch
%% OUTPUT
%%  AllPosPho:   double array, index of pho positions
patches=length(positions(:,1));
AllPosPho=[];
for i=1:1:patches
    AllPosPho=[AllPosPho,positions(i,2):1:positions(i,3)];
end

end