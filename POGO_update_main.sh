mkdir POGO_v1
cp data.mat POGO_v1
cp 16S_rRNA.fna POGO_v1
cp -r markergene POGO_v1
mv pogo_info.csv POGO_v1
mv gene_ids.csv POGO_v1

mkdir update
cd update
mkdir intermediate_archive

# 01.get genbank_list.txt from NCBI
../codes/get_genbank_list.sh

# 02.get new_genomes.list
matlab -nodisplay -nosplash -nodesktop -r "addpath ../codes; run('get_new_genome_names.m'); exit";
mv genbank_list.txt intermediate_archive

# 03.get gbk files
../codes/get_new_genome_gbks.sh

# 04.get fasta files
../codes/get_new_genome_fasta.sh

# 05.extract 16S and CDS files
../codes/get_16S_and_CDS.sh

# 06.align 16S (include parallel options)
matlab -nodisplay -nosplash -nodesktop -r "addpath ../codes; run('align_16S.m'); exit";

# 07.parse for genome pairs whose 16S are >80% similar
matlab -nodisplay -nosplash -nodesktop -r "addpath ../codes; run('get_genome_pairs.m'); exit";

# 08.get taxonomy
while read name; do
        echo $name
        uid=$(echo $name | tr "_" "\n" | tail -1)
        ../codes/extract_taxonomy.py $uid
done < new_genome_uids.list
matlab -nodisplay -nosplash -nodesktop -r "addpath ../codes; run('get_taxonomy.m'); exit";
mv taxonomy.txt intermediate_archive
mv new_genome_names.list intermediate_archive
mv new_genome_uids.list intermediate_archive

# 08.build blastdb (prot from CDS & nucl from fasta)
../codes/build_blastdb.sh

# 09.blast (include parallel options)
mkdir /scratch/yl497/1v2 /scratch/yl497/2v1
cp blastdb/* /scratch/yl497/blastdb/
cp CDS/* /scratch/yl497/CDS/
mkdir script
count=0
while read name; do
        count=$((count+1))
        uidA=$(echo $name | tr "," "\n" | head -n1)
        uidB=$(echo $name | tr "," "\n" | tail -n1)
        sed -e "s|uidA|$uidA|g" ../codes/blast_job_template.script > temp.script
        sed -e "s|uidB|$uidB|g" temp.script > script/blast$count.script
        rm temp.script
done < genome_pair.list
for i in `ls script/blast*.script`;do qsub $i;done

# 10.get orthologs
mkdir /scratch/yl497/ortholog
while read name; do
        uidA=$(echo $name | tr "," "\n" | head -n1)
        uidB=$(echo $name | tr "," "\n" | tail -n1)
        sort -u -k1,1 /scratch/yl497/1v2/${uidA}_vs_${uidB}.blasttabularout | awk '{print $1 "\t" $2}' | sort > temp1
        sort -u -k1,1 /scratch/yl497/2v1/${uidB}_vs_${uidA}.blasttabularout | awk '{print $2 "\t" $1}' | sort > temp2
        comm -12 temp1 temp2 | sed "s|s||g" > /scratch/yl497/ortholog/${uidA}_vs_${uidB}.ortholog
        rm temp1 temp2
done < genome_pair.list

# 11.realign orthologs
mkdir realign
count=0
while read name; do
        count=$((count+1))
        uidA=$(echo $name | tr "," "\n" | head -n1)
        uidB=$(echo $name | tr "," "\n" | tail -n1)
        sed -e "s|uidA|$uidA|g" ../codes/realign_job_template.script > temp
        sed -e "s|uidB|$uidB|g" temp > script/realign$count
        rm temp
done < genome_pair.list
for i in `ls script/realign*`;do qsub $i;done

# 12.read realign results and compute pairwise metrics
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('realign_summary.m'); exit";

# 13.blast ref_markergene.fasta against each genome
../codes/blast_ref_markergene.sh

# 14.parse blast output and write markergene_temp fasta files
mkdir markergene_temp
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('parse_ref_markergene.m'); exit";

# 15.blast markergene_temp fasta files against ref_markergene
../codes/blast_against_ref_markergene.sh

# 16.parse blast output and write markergene fasta files
mkdir markergene
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('parse2_ref_markergene.m'); exit";
mv markergene_temp/ intermediate_archive/

# 17.merge existing database with update database and compute fluidity
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('merge_databases.m'); exit";
cat ../16S_rRNA.fna 16S_rRNA.fna > temp; mv temp 16S_rRNA.fna
mv blastdb/* ../blastdb/; rm -r blastdb
mv gbk/* ../gbk/; rm -r gbk
mv fasta/* ../fasta/; rm -r fasta
mv CDS/* ../CDS/; rm -r CDS
mv blast_files/1v2/* ../blast_files/1v2/;
mv blast_files/2v1/* ../blast_files/2v1/;
mv blast_files/blastRef/* ../blast_files/blastRef/;
mv blast_files/blastAgainstRef/* ../blast_files/blastAgainstRef/;
rm -r blast_files
mkdir temp
for i in `ls markergene`;do cat ../markergene/$i markergene/$i >> temp/$i;done
rm -r markergene; mv temp markergene
mv genome_pair.list intermediate_archive

# 18.compute pairwise identity of markergene
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('align_markergene.m'); exit";

# 19.generate files for POGO-DB update
matlab -nodisplay -nosplash -nodesktop -nojvm -r "addpath ../codes; run('generate_POGO_files.m'); exit";

# 20.rename folder as POGO_version
cp -r markergene/ ..
cp data.mat ..
cp 16S_rRNA.fna ..
cd ..
mv update/ POGO_v2/

# 20.MySQL
mysql -u root --local-infile=1 -p data_new
# CREATE TABLE taxonomy_new LIKE taxonomy; 
# LOAD DATA LOCAL INFILE 'POGO_taxonomy.csv' INTO TABLE taxonomy_new FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
# CREATE TABLE taxonomy_api_new SELECT * FROM taxonomy_new;
# CREATE TABLE data_new LIKE data; 
# LOAD DATA LOCAL INFILE 'POGO_data.csv' INTO TABLE data_new FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
# CREATE TABLE data_and_taxonomy_new LIKE data_and_taxonomy; 
# LOAD DATA LOCAL INFILE 'POGO_data_and_taxonomy.csv' INTO TABLE data_and_taxonomy_new FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
# RENAME TABLE taxonomy TO taxonomy_bak;
# RENAME TABLE taxonomy_api TO taxonomy_api_bak;
# RENAME TABLE data TO data_bak;
# RENAME TABLE data_and_taxonomy TO data_and_taxonomy_bak;
# RENAME TABLE taxonomy_new TO taxonomy;
# RENAME TABLE taxonomy_api_new TO taxonomy_api;
# RENAME TABLE data_new TO data;
# RENAME TABLE data_and_taxonomy_new TO data_and_taxonomy;

# POGO update completed. Yeah!

