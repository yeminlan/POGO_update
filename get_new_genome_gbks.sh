mkdir gbk -p

#for each name in the list download the corresponding genbank files
while read name; do
	echo $name
	mkdir "gbk/$name"
	wget -nd -r -l 1 -A gbk "ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/$name/" -P "gbk/$name/"
        uid=$(echo $name | tr "_" "\n" | tail -1)
        cat gbk/$name/*.gbk > gbk/$uid.gbk
        rm -r gbk/$name
done < new_genome_names.list

