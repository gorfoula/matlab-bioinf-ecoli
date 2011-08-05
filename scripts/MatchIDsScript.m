
[FileNames] = ReadTable('JoinFiles/OtherStudies/NewStudies/Files.txt');
files=size(FileNames,1);

for f=1:1:files
    display(FileNames{f,1});
    param=FileNames{f,2};
    [start_idx, end_idx, extents, matches, tokens, names, splits]=regexp(param,'[ ]');
    param_d=CellTable2Double(splits);
    param_d_resh=reshape(param_d,2,length(param_d)/2)';
    FishLines('JoinFiles/uniprot-organism_83333.txt',FileNames{f,1},param_d_resh,[1 0]);
end