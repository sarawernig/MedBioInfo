#! /usr/bin/env python3

'''
Made in Python 3.6
Author: Sara Wernig Zorc
Usage: ./motif_preference_2.py
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


with open("sites/MA0106.1.sites") as handle:
     p53 = motifs.read(handle, "sites")

motif = p53.degenerate_consensus

with open("motif_result_p53.txt","w") as f:
	for seq_record in SeqIO.parse('gencode.v26.lncRNA_transcripts.fa','fasta'):
		f.write(">" + str(seq_record.id) + "\n")
		result=SeqUtils.nt_search(str(seq_record), motif)
		f.write(str(result) + "\n")

##

with open("sites/MA0001.1.sites") as handle:
     AGL3 = motifs.read(handle, "sites")

motif = AGL3.degenerate_consensus

with open("motif_result_AGL3.txt","w") as f:
	for seq_record in SeqIO.parse('gencode.v26.lncRNA_transcripts.fa','fasta'):
		f.write(">" + str(seq_record.id) + "\n")
		result=SeqUtils.nt_search(str(seq_record), motif)
		f.write(str(result) + "\n")

##

with open("sites/MA0024.1.sites") as handle:
     E2F1 = motifs.read(handle, "sites")

motif = E2F1.degenerate_consensus

with open("motif_result_E2F1.txt","w") as f:
	for seq_record in SeqIO.parse('gencode.v26.lncRNA_transcripts.fa','fasta'):
		f.write(">" + str(seq_record.id) + "\n")
		result=SeqUtils.nt_search(str(seq_record), motif)
		f.write(str(result) + "\n")



#For giving updates to the user on the progress of the program
#sys.stdeer.write()

#Gives number of occuranes of the motif
#open('windslid_10bp_out.txt', 'r').read().find('AGACATGCCT')
#TP53: 323910825

#For looping through frequency matrix's
#fn = open("./pfm_vertebrates.txt")
#for m in motifs.parse(fn, "jaspar"):

#Import matrix as a numpy array???
#Assign scores 1-4, based on numbers in the frequency matrix
