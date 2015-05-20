load data.mat

f=fopen('POGO_taxonomy.csv','wt');
fprintf(f,'id\tGenome\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tSuperkingdom\n');
for i=1:numel(genomes)
  fprintf(f,[num2str(i),'\t']);
  fprintf(f,[genomes{i},'\t']);
  fprintf(f,[taxonomy.phylum{i},'\t']);
  fprintf(f,[taxonomy.class{i},'\t']);
  fprintf(f,[taxonomy.order{i},'\t']);
  fprintf(f,[taxonomy.family{i},'\t']);
  fprintf(f,[taxonomy.genus{i},'\t']);
  fprintf(f,[taxonomy.species{i},'\t']);
  fprintf(f,[taxonomy.superkingdom{i},'\n']);
end
fclose(f);
clear f ans i

f=fopen('POGO_data.csv','wt');
fprintf(f,'id\tgenome_id1\tgenome_id2\tnumber_of_genes1\tnumber_of_genes2\tfile1v2\tfile2v1\torthologs_criterion1\torthologs_criterion2\tAverage_Amino_Acid_Identity\tGenomic_Fluidity\t16S_rRNA\tArgS\tCdsA\tCoaE\tCpsG\tDnaN\tEfp \tExo\tFfh \tFtsY\tFusA\tGlnS\tGlyA\tGroL\tHisS\tIleS\tInfA\tInfB\tKsgA\tLeuS\tMap \tMetG\tNrdA\tNusG\tPepP\tPheS\tPheT\tProS\tPyrG\tRecA\tRplA\tRplB\tRplC\tRplD\tRplE\tRplF\tRplJ\tRplK\tRplM\tRplN\tRplP\tRplR\tRplV\tRplX\tRpoA\tRpoB\tRpoC\tRpsB\tRpsC\tRpsD\tRpsE\tRpsG\tRpsH\tRpsI\tRpsJ\tRpsK\tRpsL\tRpsM\tRpsN\tRpsO\tRpsQ\tRpsS\tSecY\tSers\tThrS\tTmk \tTopA\tTrpS\tTruB\tTrxA\tTrxB\tTufB\tTyrS\tValS\n');
for i=1:size(genome_pair,1)
  fprintf(f,[num2str(i),'\t']);
  t1=find(strcmp(uid,genome_pair{i,1}));
  t2=find(strcmp(uid,genome_pair{i,2}));
  t1=t1(end);
  t2=t2(end);
  fprintf(f,[num2str(t1),'\t']);
  fprintf(f,[num2str(t2),'\t']);
  fprintf(f,[num2str(numberofgenes(t1)),'\t']);
  fprintf(f,[num2str(numberofgenes(t2)),'\t']);
  fprintf(f,['1v2/',genome_pair{i,1},'_vs_',genome_pair{i,2},'.blasttabularout.bz2\t']);
  fprintf(f,['2v1/',genome_pair{i,2},'_vs_',genome_pair{i,1},'.blasttabularout.bz2\t']);
  fprintf(f,[num2str(orthologs1(i)),'\t']);
  fprintf(f,[num2str(orthologs2(i)),'\t']);
  fprintf(f,[num2str(AAI(i),'%.8f'),'\t']);
  fprintf(f,[num2str(fluidity(i),'%.8f')]);
  for j=1:numel(gene_name)
    fprintf(f,['\t',num2str(gene_identity(i,j),'%.8f')]);
  end
  fprintf(f,'\n');
end
fclose(f);
clear f ans i j t1 t2

f=fopen('POGO_data_and_taxonomy.csv','wt');
fprintf(f,'id\tgenome_id1\tgenome_id2\tnumber_of_genes1\tnumber_of_genes2\tfile1v2\tfile2v1\torthologs_criterion1\torthologs_criterion2\tAverage_Amino_Acid_Identity\tGenomic_Fluidity\t16S_rRNA\tArgS\tCdsA\tCoaE\tCpsG\tDnaN\tEfp \tExo\tFfh \tFtsY\tFusA\tGlnS\tGlyA\tGroL\tHisS\tIleS\tInfA\tInfB\tKsgA\tLeuS\tMap \tMetG\tNrdA\tNusG\tPepP\tPheS\tPheT\tProS\tPyrG\tRecA\tRplA\tRplB\tRplC\tRplD\tRplE\tRplF\tRplJ\tRplK\tRplM\tRplN\tRplP\tRplR\tRplV\tRplX\tRpoA\tRpoB\tRpoC\tRpsB\tRpsC\tRpsD\tRpsE\tRpsG\tRpsH\tRpsI\tRpsJ\tRpsK\tRpsL\tRpsM\tRpsN\tRpsO\tRpsQ\tRpsS\tSecY\tSers\tThrS\tTmk \tTopA\tTrpS\tTruB\tTrxA\tTrxB\tTufB\tTyrS\tValS\tgenome1_name\tgenome1_Phylum\tgenome1_class\tgenome1_order\tgenome1_family\tgenome1_genus\tgenome1_species\tgenome1_superkingdom\tgenome2_name\tgenome2_phylum\tgenome2_class\tgenome2_order\tgenome2_family\tgenome2_genus\tgenome2_species\tgenome2_superkingdom\n');
for i=1:size(genome_pair,1)
  fprintf(f,[num2str(i),'\t']);
  t1=find(strcmp(uid,genome_pair{i,1}));
  t2=find(strcmp(uid,genome_pair{i,2}));
  t1=t1(end);
  t2=t2(end);
  fprintf(f,[num2str(t1),'\t']);
  fprintf(f,[num2str(t2),'\t']);
  fprintf(f,[num2str(numberofgenes(t1)),'\t']);
  fprintf(f,[num2str(numberofgenes(t2)),'\t']);
  fprintf(f,['1v2/',genome_pair{i,1},'_vs_',genome_pair{i,2},'.blasttabularout.bz2\t']);
  fprintf(f,['2v1/',genome_pair{i,2},'_vs_',genome_pair{i,1},'.blasttabularout.bz2\t']);
  fprintf(f,[num2str(orthologs1(i)),'\t']);
  fprintf(f,[num2str(orthologs2(i)),'\t']);
  fprintf(f,[num2str(AAI(i),'%.8f'),'\t']);
  fprintf(f,[num2str(fluidity(i),'%.8f'),'\t']);
  for j=1:numel(gene_name)
    fprintf(f,[num2str(gene_identity(i,j),'%.8f'),'\t']);
  end
  fprintf(f,[genomes{t1},'\t']);
  fprintf(f,[taxonomy.phylum{t1},'\t']);
  fprintf(f,[taxonomy.class{t1},'\t']);
  fprintf(f,[taxonomy.order{t1},'\t']);
  fprintf(f,[taxonomy.family{t1},'\t']);
  fprintf(f,[taxonomy.genus{t1},'\t']);
  fprintf(f,[taxonomy.species{t1},'\t']);
  fprintf(f,[taxonomy.superkingdom{t1},'\t']);
  fprintf(f,[genomes{t2},'\t']);
  fprintf(f,[taxonomy.phylum{t2},'\t']);
  fprintf(f,[taxonomy.class{t2},'\t']);
  fprintf(f,[taxonomy.order{t2},'\t']);
  fprintf(f,[taxonomy.family{t2},'\t']);
  fprintf(f,[taxonomy.genus{t2},'\t']);
  fprintf(f,[taxonomy.species{t2},'\t']);
  fprintf(f,[taxonomy.superkingdom{t2},'\n']);
end
fclose(f);
clear f ans i j t1 t2

