load data.mat new_uid
load ../data.mat gene_name
gene_name = gene_name(2:end);%excluding 16S

for i=1:numel(new_uid)
    %read blastAgainstRef output
    gene_exist=zeros(numel(gene_name),1);
    blast=myblastreadlocal(['blast_files/blastAgainstRef/',new_uid{i},'.blasttabularout'],8);
    for j=1:numel(blast)
        if strcmp(blast(j).Query,blast(j).Hits(1).Name)
            gene_exist(strcmp(gene_name,blast(j).Query))=1;
        end
    end
    clear blast j
    gene_exist=gene_name(logical(gene_exist));

    %write markergene fasta files in markergene/
    f1=myfastaread(['markergene_temp/',new_uid{i},'.fasta']);
    f1=f1(ismember({f1.Header},gene_exist));
    for j=1:numel(f1)
        f=fopen(['markergene/',f1(j).Header,'.fna'],'at');
        fprintf(f,['>',new_uid{i},'\n']);
        fprintf(f,[f1(j).Sequence,'\n']);
        fclose(f);
    end
    clear j f1 f ans gene_exist
end
%clear i
