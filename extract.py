#!/usr/bin/python
from Bio import SeqIO
import sys
import os.path

genome = SeqIO.parse('gbk/' + sys.argv[1] + '.gbk', "genbank")
cds = []
s16 = []

if len(sys.argv) < 5:
	product_16s = '16s ribosomal rna'
else:
	product_16s = sys.argv[4]

product_16s = product_16s.lower() 

print "Searching 16s tag:", product_16s
for record in genome:
	for feature in record.features:
		
		# extract 16s
		if 'product' in feature.qualifiers:
			product=feature.qualifiers['product'][0].lower()
			if product == product_16s:
                                if len(feature.extract(record).seq) >= 1000:
                                        if len(feature.extract(record).seq) <= 1800:
				                s16.append(feature.extract(record).seq)

		# extract all CDS protiens
		if 'translation' in feature.qualifiers:
			cds.append(feature.qualifiers['translation'][0])

# get basename with no extension
fn,_ = os.path.splitext(sys.argv[1])

s16 = set(s16)

with open('16S/' + sys.argv[1] + '.fna', 'w') as fh:
	if len(s16) != 0:
		for i, it in enumerate(s16):
                        fh.write(">")
                        fh.write(os.path.basename(fn))
                        fh.write('_' + str(i+1))
                        fh.write('\n')
			fh.write(str(it))
			fh.write('\n');

with open('CDS/' + sys.argv[1] + '.faa', 'w') as fh:
	if len(cds) != 0:
		for i, it in enumerate(cds):
			fh.write('>s' + str(i+1))
			fh.write('\n')
			fh.write(str(it))
			fh.write('\n')
