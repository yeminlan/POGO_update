mkdir 16S CDS

while read name; do
        echo $name
        uid=$(echo $name | tr "_" "\n" | tail -1)
        ../codes/extract.py $uid
done < new_genome_names.list

cat 16S/*.fna > 16S_rRNA.fna
rm -r 16S


