#!/usr/bin/python
from Bio import SeqIO
import sys
import os.path

genome = next(SeqIO.parse('gbk/' + sys.argv[1] + '.gbk', "genbank"))
s = genome.annotations['taxonomy']
s.append(genome.annotations['organism'])

with open('taxonomy.txt', 'a+') as fh:
	fh.write(sys.argv[1] + "\t")
	for i in s:
		fh.write(i + "\t")
	fh.write('\n');
