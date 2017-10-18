#! /usr/bin/env python3
#-*- coding: utf-8 -*-

'''
Made in Python 3.6
Author: Sara Wernig Zorc
Copyright (C) 2017 Sara Wernig Zorc
Usage: ./motif_preference_2.py input_file.txt outputfile.txt
'''

'''
import sys		#For standard error report and system arguments
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
all_positions=[]
exon_positions = {}
gene_positions = {}
motif_positions = {}

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

'''Module find_lincRNA, gives all lincRNA exons in the GTF file, with transcriptID,
Chromosome, star and end positions of an exon in a dicitonary'''

def find_lincRNA_exon(annotation_file):
	with open(annotation_file, "r") as gtf:
		for line in gtf:
			if line.startswith('#'):
				continue

			else:
				b=line.split(sep='\t')
				c=b[-1].split(sep=';')
				d=c[1].split(sep=' ')

				#print(b)
				if 'exon' in b[2]:
					#print(line)
					trans_id = str(d[-1]).strip("\" ")
					Chromosome=b[0]
					start=int(b[3])
					stop=int(b[4])
					exon_name=c[6].strip()
					one_exon = (exon_name,Chromosome,start,stop)
					list_exons = [one_exon]
					#print(one_exon)

					if trans_id in exon_positions.keys():
						exon_positions[trans_id].append(one_exon)
					else:
						empty_values = exon_positions.get(trans_id,None)
						exon_positions[trans_id] = list_exons

	#print(exon_positions)
	return exon_positions

def find_lincRNA_transcript(annotation_file):
	with open(annotation_file, "r") as gtf:
		for line in gtf:
			if line.startswith('#'):
				continue

			else:
				b=line.split(sep='\t')
				c=b[-1].split(sep=';')
				d=c[1].split(sep=' ')

				if 'transcript' in b[2]:
					#print(line)
					trans_id = str(d[-1]).strip("\" ")
					Chromosome=b[0]
					start=int(b[3])
					stop=int(b[4])
					list_gene = [Chromosome,start,stop]
					#print(list_gene)

					if trans_id in gene_positions.keys():
						gene_positions[trans_id].append(list_gene)
					else:
						empty_values = gene_positions.get(trans_id,None)
						gene_positions[trans_id] = list_gene

	#print(gene_positions)
	return gene_positions

def lincRNA_sequence(input_fasta,annotation_file,):
	for line in find_lincRNA(annotation_file):
		a=str(genome_sequence(input_fasta))
		start=int(line[1])
		stop=int(line[2])
		lincRNA_seq = a[start:stop]
		pro_start = (start - 2500)
		pro_stop = (start - 1)
		promoter_seq = a[pro_start:pro_stop]
		key = line
		#print(lincRNA_seq)
		whole_sequence = dic.get(key,None,None)
		dic[key] = (lincRNA_seq, promoter_seq)
		print(dic)
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

'''Threshold of log-odds 7 = 100x more likely to occur in motif
than random background.
Negative positions are on - strand, positive positions are on + strand.
A highly selective motif should only match once (or zero times)
in each sequence tested.'''

def motif_search(input_jaspar,lincRNA_seq,output_motif):
	with open (input_jaspar,"r") as fm, open(lincRNA_seq,"r") as lincs, open(output_motif,"w") as f:
		for m in motifs.parse(fm,"jaspar"):
			for lincseq in SeqIO.parse(lincs,'fasta',alphabet=IUPAC.unambiguous_dna):
				for pos, score in m.pssm.search(lincseq.seq, threshold=7.0):
					d=lincseq.id.split("|")
					linc_id = d[1].strip("\" ")
					position = abs(pos)
					final_position = (int(position)/int(d[-2]))
					all_positions.append(final_position)
					#f.write("For lincRNA" + "\t" + str(d[0]) + "\t" + str(d[1]) + "\t" + "total length:" + "\t" + str(d[-2]) + "\n")
					#f.write("Motif" + "\t" + str(m.name) + "\t" + "binds at position" + "\t" + str(pos) + "\t" + " with score: " + "\t" + str(score) + "\n")

	return all_positions

def calculation():
	import statistics
	mode=statistics.mode(motif_search())
	mean=statistics.mean(motif_search())
	stdev=statistics.stdev(motif_search())
	var=statistics.variance(motif_search())
	print("Mode:\t" + mode + "\t" + "Mean:\t"+ mean + "Standard diviation::\t" + stdev + "Variance:\t", var)

def motif_search_2(input_jaspar,lincRNA_seq):
	with open (input_jaspar,"r") as fm, open(lincRNA_seq,"r") as lincs:
		for m in motifs.parse(fm,"jaspar"):
			for lincseq in SeqIO.parse(lincs,'fasta',alphabet=IUPAC.unambiguous_dna):
				d = lincseq.id.split("|")
				linc_id = d[0].strip("\" ")
				#print(linc_id)
				for pos, score in m.pssm.search(lincseq.seq):
					one_position = int(abs(pos))
					list_positions = [one_position]
					#print(pos)
					if linc_id in motif_positions.keys():
						motif_positions[linc_id].append(one_position)
					else:
						empty_values = motif_positions.get(linc_id,None)
						motif_positions[linc_id] = list_positions

	#print(motif_positions)
	return motif_positions

def main_preference_of_binding(input_jaspar,lincRNA_seq,annotation_file):
	a = motif_search_2(input_jaspar,lincRNA_seq)
	b = find_lincRNA_transcript(annotation_file)
	c = find_lincRNA_exon(annotation_file)

	for key, value in a.items():
		print(key,value)

	for key, value in b.items():
		print(key,value)

	for key, value in c.items():
		print(key,value)



if __name__ == '__main__':
	#genome_sequence("input/hg38.fa")
	#find_lincRNA_transcript(annotation_file="input/test_file.gtf")
	#find_lincRNA_gene(annotation_file="input/test_file.gtf")
	#lincRNA_sequence(input_fasta="input/hg38.fa",annotation_file="input/gencode.v26.annotation.gtf")
	#make_consensus("input/frequency_matrixes_vertebrates_nr_JASPAR.txt","input/degConsensusSeq.txt")
	#calculation()
	#motif_search("input/JASPAR2018.txt","input/gencode.v27.lncRNA_transcripts.fa","output/test_motif_result.positions.txt")
	#motif_search_2("input/JASPAR2018.txt","input/test.fa")
	main_preference_of_binding("input/JASPAR2018.txt","input/test.fa","input/test_file.gtf")

#For giving updates to the user on the progress of the program
#sys.stdeer.write()
