#! /usr/bin/env python3

'''
Made in Python 3.6
Author: Sara Wernig Zorc
Copyright (C) 2017 Sara Wernig Zorc
Usage: ./make_consensus_from_matrix.py input_file.txt outputfile.txt
'''

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqUtils

def make_consensus(input_matrix,output_consensus):
	with open(output_consensus,"a") as f, open(input_matrix,"r") as matrix:
		for m in motifs.parse(matrix,"jaspar"):
			f.write(">" + str(m.base_id) + "\n")
			f.write(str(m.degenerate_consensus) + "\n")


if __name__ == '__main__':
	make_consensus("input/frequency_matrixes_vertebrates_nr_JASPAR.txt","input/degenerateConsensusSequence.txt")
