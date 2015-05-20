load data.mat
orthologs1 = zeros(size(genome_pair,1),1);
orthologs2 = zeros(size(genome_pair,1),1);
AAI = zeros(size(genome_pair,1),1);
error=[];
for i=1:size(genome_pair,1)
    try
        eval(['load realign/',genome_pair{i,1},'_vs_',genome_pair{i,2},'.mat hits'])
        if ~isempty(hits)
            filter1=((hits(:,3)>0.3)&(hits(:,4)>0.7));
            filter2=((hits(:,3)>0.1)&(hits(:,4)>0.5));
            orthologs1(i) = sum(filter1);
            orthologs2(i) = sum(filter2);
            AAI(i) = mean(hits(filter1,3));
            clear filter1 filter2
        else
            error=[error;i];
        end
    catch
        error=[error;i];
    end
    if mod(i,1000)==1
        save data2.mat
    end
end
clear i
save data2.mat

%rerun all genome pairs in error and make sure there's no errors other than 0-hit cases, then:
clear error
new_genome_pair = genome_pair;
new_numberofgenes = numberofgenes;
new_rrna_identity = rrna_identity;
new_AAI = AAI;
new_orthologs1 = orthologs1;
new_orthologs2 = orthologs2;
clear genome_pair numberofgenes rrna_identity AAI orthologs1 orthologs2
save data.mat %replace the file

