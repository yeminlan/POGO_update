load data.mat new_uid

for i=1:numel(new_uid)
    %read blastRef output
    blast=myblastreadlocal(['blast_files/blastRef/',new_uid{i},'.blasttabularout'],8);
    genome=myfastaread(['fasta/',new_uid{i},'.fna']);
    for j=1:numel(genome)
        t=regexp(genome(j).Header,' ','split');
        genome_name{j}=t{1};
    end
    clear j t
    genome_seq={genome.Sequence};
    clear genome

    %write markergene_temp fasta files
    f=fopen(['markergene_temp/',new_uid{i},'.fasta'],'wt');
    for j=1:numel(blast)
        t2=blast(j).Hits(1).HSPs(1).SubjectIndices;
        if t2(1)<=t2(2)
            t1=genome_seq{strcmp(genome_name,blast(j).Hits(1).Name)}(t2(1):t2(2));
        else
            t1=myseqrcomplement(genome_seq{strcmp(genome_name,blast(j).Hits(1).Name)}(t2(2):t2(1)));
        end
        fprintf(f,['>',blast(j).Query,'\n']);
        fprintf(f,[t1,'\n']);
        clear t1 t2 t3
    end
    fclose(f);
    clear j genome_name genome_seq
end
clear i


