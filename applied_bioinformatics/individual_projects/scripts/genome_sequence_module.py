#! /usr/bin/env python3
#-*- coding: utf-8 -*-

'''
Made in Python 3.6
Author: Sara Wernig Zorc
Copyright (C) 2017 Sara Wernig Zorc
Usage: ./motif_preference_2.py input_file.txt outputfile.txt
'''

#Biopython modules
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqUtils

genome=[]

def genome_sequence(input_fasta):
	with open(input_fasta,'r') as fa:
		for line in fa:
			if line.startswith('>'):
				genome.append(line[:1]+'Chr'+line[1:2]+'\n')
				#print(genome)
			else:
				newline=line.strip()
				genome.append(newline)
				#print(genome)

	genome_seq = ''.join(genome)
	#print(genome_seq)

if __name__ == '__main__':
	genome_sequence("input/hg38.fa")
