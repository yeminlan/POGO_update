mkdir blast_files/blastRef
while read uid; do
        tblastn -query ../codes/ref_markergene.fasta -db blastdb/$uid -out blast_files/blastRef/$uid.blasttabularout -outfmt 6
done < intermediate_archive/new_genome_uids.list

