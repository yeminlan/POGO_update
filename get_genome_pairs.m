load data.mat
X=[new_genomes,new_uid];%going to update this
clear new_genomes new_uid

load intermediate_archive/16S_similarity.mat new_uid all_uid similarity
for i = 1:numel(new_uid)
    new_genomes(i,1) = X(strcmp(X(:,2),new_uid{i}),1);
end
clear i X

[ind1,ind2] = find(similarity>=0.8);
genome_pair = [new_uid(ind1),all_uid(ind2)];
rrna_identity = similarity(find(similarity>=0.8));
clear all_uid ind1 ind2 similarity

%remove self comparisons
ind=strcmp(genome_pair(:,1),genome_pair(:,2));
genome_pair(ind,:)=[];
rrna_identity(ind,:)=[];
clear ind

% write to file genome_pair.list
f=fopen('genome_pair.list','wt');
for i = 1:size(genome_pair,1)
    fprintf(f,[genome_pair{i,1},',',genome_pair{i,2},'\n']);
end
fclose(f);
clear f ans i

% get numberofgenes for each new_genome
numberofgenes = zeros(numel(new_uid),1);
for i = 1:numel(new_uid)
    numberofgenes(i,1) = numel(myfastaread(['CDS/',new_uid{i},'.faa']));
end
clear i

save data.mat %replace the orginal one
