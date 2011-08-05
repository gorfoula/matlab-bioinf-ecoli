TheoReadTable=ReadTable('Data/DataSetGO/DB_SeptV4.txt');

proteins=size(TheoReadTable,1);
OF_F='Data/DataSetGO/DB_SeptV5.txt';
links=cell(proteins,1);

links(1)=TheoReadTable(1,9);
for i=2:proteins
    links{i}=['=HYPERLINK("http://www.uniprot.org/uniprot/?query="&','B',num2str(i),'&"+ECOLI&sort=score";B',num2str(i),')'];
end
FileWriteTable(OF_F,[TheoReadTable(:,1:8) links TheoReadTable(:,10:end)],[],'w');