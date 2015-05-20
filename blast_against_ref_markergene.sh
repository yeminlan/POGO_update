mkdir blast_files/blastAgainstRef
while read uid; do
        tblastx -query markergene_temp/$uid.fasta -db ../codes/ref_markergene -out blast_files/blastAgainstRef/$uid.blasttabularout -outfmt 6
done < intermediate_archive/new_genome_uids.list

