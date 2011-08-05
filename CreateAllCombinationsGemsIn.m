function []=CreateAllCombinationsGemsIn()

Categories={'Cytoplasmic' 'Secreted' 'Periplasmic' 'b_barrel' 'Lipoproteins'};
NumCatg=length(Categories);

dirTrain='Dataset/TRAIN/';
dirTest='Dataset/TEST/';
outSubFolder='PercTest/';
fastaSubFolder='temp/';

for i=1:1:NumCatg
    for j=i+1:1:NumCatg
        
        filename_SP=['SP',Categories{i},'VS',Categories{j}];
        filename_MAT=[Categories{i},'VS',Categories{j}];
        filename_GROUP=['SP',Categories{i},'VS',Categories{j},'_20119'];
              
        s=fopen('Input_files.txt','w');
        fprintf(s,'%s\n%s',[dirTrain,Categories{i},'_mn.txt'],[dirTrain,Categories{j},'_mn.txt']); fclose all;
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_SP,'.gemsin'],'bin',100,-1,0,0,1,[fastaSubFolder,filename_SP,'.fasta'],0);
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_MAT,'.gemsin'],'bin',100,1,0,0,1,[fastaSubFolder,filename_MAT,'.fasta'],0);
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_GROUP,'.gemsin'],'bin',100,-1,0,0,1,[fastaSubFolder,filename_GROUP,'.fasta'],20119);

        s=fopen('Input_files.txt','w');
        fprintf(s,'%s\n%s',[dirTest,Categories{i},'_mn_test.txt'],[dirTest,Categories{j},'_mn_test.txt']);fclose all;
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_SP,'.test'],'bin',100,-1,0,0,0,'temp.fasta',0);
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_MAT,'.test'],'bin',100,1,0,0,0,'temp.fasta',0);
        AllCategoriesCoding('Input_files.txt',[outSubFolder,filename_GROUP,'.test'],'bin',100,-1,0,0,0,'temp.fasta',20119);
        
    end

end


end