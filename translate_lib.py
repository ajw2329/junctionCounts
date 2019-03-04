import splice_lib
import sys

def translate_seq(seq):
	'''
		Code borrowed from http://www.petercollingridge.co.uk/python-bioinformatics-tools/codon-table
	'''

	bases = ['t', 'c', 'a', 'g']

	codons = [a+b+c for a in bases for b in bases for c in bases]

	amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'

	codon_table = dict(zip(codons, amino_acids))

	seq = seq.lower().replace('\n', '').replace(' ', '')

	peptide = ''
  
	for i in xrange(0, len(seq), 3):

		codon = seq[i: i+3]
		amino_acid = codon_table.get(codon, '*')
		if amino_acid != '*':

			peptide += amino_acid

		else:

			break
                 
	return peptide

