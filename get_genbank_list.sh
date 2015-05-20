#!/bin/bash
curl -l ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/ --user anonymous:yl497@drexel.edu | grep "_uid" > genbank_list.txt

