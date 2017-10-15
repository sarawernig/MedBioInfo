#! /usr/bin/env python3
#-*- coding: utf-8 -*-

'''
Made in Python 3.6
Author: Sara Wernig Zorc
Copyright (C) 2017 Sara Wernig Zorc
Usage: ./motif_preference_2.py input_file.txt outputfile.txt
'''

'''
import sys		#For standard error report to the user and debugging
import numpy		#Import module for matrixes
import re		#Import module for regular expressions
'''

#Biopython modules
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqUtils

genome=[]
lincRNAs=[]
dic={}

'''Module genome_sequence, gives the whole genome sequence as one line,
seperated by >Chromosome name'''

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
	return genome_seq

def find_lincRNA(annotation_file):
	with open(annotation_file) as gtf:
		for line in gtf:
			if line.startswith('#'):
				continue
			else:
				b=line.split(sep='\t')
				c=b[-1].split(sep=';')
				if '\ttranscript\t' in line and 'lincRNA' in c[1]:
					#print(line)
					d=c[0].split(sep=' ')
					trans_id = d[2]
					lincRNA = (b[0],b[3],b[4],trans_id)
					#print(lincRNA)
					lincRNAs.append(lincRNA)
	return lincRNAs

def lincRNA_sequence(find_lincRNA):
	for line in lincRNAs:
		linRNA_seq = genome_seq[line[1]:line[2]]
		#print(linRNA_sequence)
		linRNA_seq=dic.get(lincRNAs,('','',''))
		#print(dic)
	return dic

def make_consensus(input_matrix,output_consensus):
	with open(output_consensus,"a") as f, open(input_matrix,"r") as matrix:
		for m in motifs.parse(matrix,"jaspar"):
			f.write(">" + str(m.base_id) + "\n")
			f.write(str(m.degenerate_consensus) + "\n")

def motif_search(input_consensus,lincRNA_seq,output_motif):
	with open(input_consensus,"r") as consensus, open(output_motif,"a") as out_file, open(lincRNA_seq,"r"):
		for line in consensus:
			if line.startswith('>'):
				motif_id = line
			else:
				motif_sequence = line

		for line in lincRNA_seq:
			if line.startswith('>'):
				b = a.split(sep='|')
				lincRNA_id = b[0]
			else:
				lincRNA_sequence = line
		print(motif_id,motif_sequence,lincRNA_id,lincRNA_sequence)
		#out_file.write(SeqUtils.nt_search(lincRNA_sequence, motif_sequence))
		#return result

def main():
	genome_sequence()
	find_lincRNA()
	lincRNA_sequence()
	make_consensus()
	motif_search()


if __name__ == '__main__':
	#genome_sequence("input/chr1.fa")
	#find_lincRNA(annotation_file="input/test_file.gtf",genome_seq)
	#lincRNA_seq(lincRNAs,genome_seq)
	#make_consensus("input/frequency_matrixes_vertebrates_nr_JASPAR.txt","input/degenerateConsensusSequence.txt")
	motif_search(input_consensus="input/degenerateConsensusSequence.txt",lincRNA_seq="input/gencode.v26.lncRNA_transcripts.fa",output_motif="output/motif_result.txt",)


#For giving updates to the user on the progress of the program
#sys.stdeer.write()
