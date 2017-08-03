#! /usr/bin/env python3

'''import sys		#For standard error report to the user and debugging
import numpy		#Import module for matrixes
import re'''		#Import module for regular expressions

#Biopython modules
from Bio import SeqIO	
from Bio.Alphabet import IUPAC		
from Bio.Seq import Seq
from Bio import motifs			
from Bio import SeqUtils

motif = Seq("ATAN", IUPAC.ambiguous_dna)			#TESTING
TP53_seq = Seq("AGACATGCCT", IUPAC.unambiguous_dna)	#TP53 10bp binding motif
CTCF_seq = Seq("CCGCGNGGNGGCAG", IUPAC.ambiguous_dna)	#CTCF 14bp binding motif

with open("motif_result.txt","w") as f:
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
