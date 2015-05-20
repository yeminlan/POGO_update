load data.mat
load ../data.mat

genome_pair = [genome_pair;new_genome_pair];
genomes = [genomes;new_genomes];
numberofgenes = [numberofgenes;new_numberofgenes];
uid = [uid;new_uid];
gene_identity = [gene_identity;new_rrna_identity,zeros(size(new_rrna_identity,1),73)];
AAI = [AAI;new_AAI];
orthologs1 = [orthologs1;new_orthologs1];
orthologs2 = [orthologs2;new_orthologs2];
clear new_genome_pair new_genomes new_numberofgenes new_uid new_rrna_identity new_AAI new_orthologs1 new_orthologs2

%compute fluidity
for i=numel(fluidity)+1:numel(AAI)
    if orthologs2(i) >= 200
        t1=numberofgenes(strcmp(uid,genome_pair{i,1}));
        t2=numberofgenes(strcmp(uid,genome_pair{i,2}));
        fluidity(i) = 1 - 2*orthologs2(i)/(t1(end)+t2(end));
        clear t1 t2
    else
        fluidity(i) = NaN;
    end
end
clear i

save data.mat %replace
