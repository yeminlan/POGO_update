mkdir fasta -p

#for each name in the list download the corresponding genbank files
while read name; do
	echo $name
	mkdir "fasta/$name"
	wget -nd -r -l 1 -A fna "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/$name/" -P "fasta/$name/"
        uid=$(echo $name | tr "_" "\n" | tail -1)
        cat fasta/$name/*.fna > fasta/$uid.fna
        rm -r fasta/$name
done < new_genome_names.list

