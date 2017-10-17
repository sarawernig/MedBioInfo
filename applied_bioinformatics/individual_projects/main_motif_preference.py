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

genome=[]
lincRNAs=[]
dic={}
lincRNA_seq=''

'''Module genome_sequence, gives the whole genome sequence as one line,
seperated by >Chromosome name'''

def genome_sequence(input_fasta):
	with open(input_fasta,'r') as fa:
		for line in fa:
			if line.startswith('>'):
				continue
				#genome.append(line[:1]+'Chr'+line[1:2]+'\n')
			else:
				newline=line.strip()
				genome.append(newline)
				#print(genome)

	genome_seq = ''.join(genome)
	#print(genome_seq[:500])
	return genome_seq

'''Module find_lincRNA, gives all lincRNAs in the GTF file, with transcriptID,
Chromosome, star and end positions in a tuple'''

def find_lincRNA(annotation_file):
	with open(annotation_file) as gtf:
		for line in gtf:
			if line.startswith('#'):
				continue
			elif '\ttranscript\t' in line:
				b=line.split(sep='\t')
				c=b[-1].split(sep=';')
				#print(line)
				if 'lincRNA' in c[2]:
					d=c[1].split(sep=' ')
					Chromosome=b[0]
					start=b[3]
					stop=b[4]
					trans_id = d[-1]
					lincRNA = (Chromosome,start,stop,trans_id)
					lincRNAs.append(lincRNA)
					#print(lincRNAs)
					return lincRNAs

def lincRNA_sequence(input_fasta,annotation_file,):
	for line in find_lincRNA(annotation_file):
		a=str(genome_sequence(input_fasta))
		start=int(line[1])
		stop=int(line[2])
		lincRNA_seq = a[start:stop]
		key = line
		#print(lincRNA_seq)
		whole_sequence = dic.get(key,None)
		dic[key] = (lincRNA_seq)
		#print(dic)
	return dic


def make_consensus(input_matrix,output_consensus):
	with open(output_consensus,"a") as f, open(input_matrix,"r") as matrix:
		for m in motifs.parse(matrix,"jaspar"):
			f.write(">" + str(m.base_id) + "\n")
			f.write(str(m.consensus) + "\n")

def make_instances(input_matrix,output_instances):
	with open(output_instances,"a") as f, open(input_matrix,"r") as matrix:
		for m in motifs.parse(matrix,"jaspar"):
			f.write(">" + str(m.base_id) + "\n")
			f.write(str(m.instances) + "\n")

def motif_search(input_jaspar,lincRNA_seq,output_motif):
	with open (input_jaspar,"r") as fm, open(lincRNA_seq,"r") as lincs, open(output_motif,"w") as f:
		for m in motifs.parse(fm,"jaspar"):
			for lincseq in SeqIO.parse(lincs,'fasta',alphabet=IUPAC.unambiguous_dna):
				for pos, score in m.pssm.search(lincseq.seq, threshold=3.0):
					d=lincseq.id.split("|")
					#print(pos,score,d[0],m.base_id)
					f.write("For lincRNA" + "\t" + str(d[0]) + "\t" + str(d[1]) + "\t" + "total length:" + "\t" + str(d[-2]) + "\n")
					f.write("Motif" + "\t" + str(m.name) + "\t" + "binds at position" + "\t" + str(pos) + "\t" + " with score: " + "\t" + str(score) + "\n")

def main():
	motif_search()


if __name__ == '__main__':
	#genome_sequence("input/chr1.fa")
	#find_lincRNA(annotation_file="input/test_file.gtf")
	lincRNA_sequence(input_fasta="input/chr1.fa",annotation_file="input/test_file.gtf")
	#make_consensus("input/frequency_matrixes_vertebrates_nr_JASPAR.txt","input/degenerateConsensusSequence.txt")
	#motif_search("input/JASPAR2018.txt","input/test.fa","output/motif_result.positions.txt")


#For giving updates to the user on the progress of the program
#sys.stdeer.write()
