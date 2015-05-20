% 01.get rna/uid

load ../data.mat uid

new_rna=fastaread('16S_rRNA.fna');
rna=fastaread('../16S_rRNA.fna');

for i=1:numel(new_rna)
    t=strfind(new_rna(i).Header,'_');
    new_rna(i).Header=new_rna(i).Header(1:t-1);
end
clear i t
for i=1:numel(rna)
    t=strfind(rna(i).Header,'_');
    rna(i).Header=rna(i).Header(1:t-1);
end
clear i t

new_uid=unique({new_rna.Header}');

all_uid=[uid;new_uid];
all_rna=[rna;new_rna];

% 02.write the file new_genome_uids.list
%It is updated so that new genomes without 16S are not in this list
f=fopen('new_genome_uids.list','wt');
for i=1:numel(new_uid)
    fprintf(f,[new_uid{i},'\n']);
end
fclose(f);
clear i ans f

% 03.align 16S, I'd break the loop into parallel jobs
similarity=zeros(numel(new_uid),numel(all_uid));
for i=1:numel(new_uid)
    parallel_align(i);%save in separate part_$i.mat files
end
clear i

% 04.put together the parallel jobs in step 03
for i=1:numel(new_uid)
    eval(['load part',num2str(i,'%.3d'),'.mat similarity'])
    Similarity(i,:)=similarity(i,:);
end
similarity=Similarity;
clear i Similarity
clear rna new_rna all_rna 

% 05.remove old_uids that are overlapping with new_uids (old genomes that have been updated)
old_uid_to_remove = find(ismember(uid,new_uid));
all_uid(old_uid_to_remove)=[];
similarity(:,old_uid_to_remove)=[];

save intermediate_archive/16S_similarity.mat


