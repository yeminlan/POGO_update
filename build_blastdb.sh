mkdir blastdb -p

while read uid; do
#        makeblastdb -in CDS/$uid.faa -dbtype prot -out blastdb/$uid -parse_seqids
        makeblastdb -in fasta/$uid.fna -dbtype nucl -out blastdb/$uid -parse_seqids
done < new_genome_uids.list
