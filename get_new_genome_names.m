%read new genome list
f=fopen('genbank_list.txt');
t_list={};
while(1)
    line=fgets(f);
    if all(line==-1)
        break
    end
    t_list=[t_list;line(1:end-1)];
end
fclose(f);
clear f ans line

%load current genome list
load ../data.mat genomes

%get new_genomes
new_genomes=t_list(~ismember(t_list,genomes));
clear t_list genomes

%get new_uid
for i=1:numel(new_genomes)
    t1=regexp(new_genomes{i},'_','split');
    new_uid{i,1}=t1{end};
end
clear t1 i

%write to new_genome_names.list
f=fopen('new_genome_names.list','wt');
for i=1:numel(new_genomes)
    fprintf(f,[new_genomes{i},'\n']);
end
fclose(f);
clear i f ans

save data.mat


