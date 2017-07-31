#! /usr/bin/env python3

import sys		#For standard error report to the user and debugging
import numpy		#Import module for matrixes
import re		#Import module for regular expressions
from Bio import SeqIO	#Import module for window sliding
from Bio.Alphabet import IUPAC		#Import module for ambigious sequences
from Bio.Seq import Seq
from Bio import motifs			#Motif search module

#Change this to a string not an output file!
#Add the position of the sequence at the ned of the header!!!

#with open("windslid_10bp_out.txt","w") as f:
#	for seq_record in SeqIO.parse("gencode.v26.lncRNA_transcripts.fa", "fasta"):
#		for i in range(len(seq_record.seq) - 9) :		#Length of sequence minus 9
#			f.write(">" + str(seq_record.id) + "\n")		#Write the name of the sequence
#			f.write(str(seq_record.seq[i:i+10]) + "\n") 	#Window sliding for 10bp motifs, moving by 1bp	

#with open("windslid_14bp_out.txt","w") as f:
#	for seq_record in SeqIO.parse("gencode.v26.lncRNA_transcripts.fa", "fasta"):
#		for i in range(len(seq_record.seq) - 13) :		#Length of sequence minus 13
#			f.write(">" + str(seq_record.id) + "\n")		#Write the name of the sequence
#			f.write(str(seq_record.seq[i:i+14]) + "\n") 	#Window sliding for 10bp motifs, moving by 1bp	



#Make a dictionary with TFs and their binding motifs??
TP53_seq = Seq("AGACATGCCT", IUPAC.ambiguous_dna)	#TP53 10bp binding motif
CTCF_seq = Seq("CCGCGNGGNGGCAG", IUPAC.ambiguous_dna)	#CTCF 14bp binding motif


#Find 10bp motifs in lncRNA sequences (only yes or no; doesn't give the position nor number of occurences)

fname='windslid_14bp_out.txt'

def check(fname, txt):
    with open(fname) as dataf:
        return any(txt in line for line in dataf)

if check('windslid_14bp_out.txt', 'CCGCGAGGAGGCAG'):
    print('true')
else:
    print('false')

#Gives number of occuranes of the motif
open('windslid_10bp_out.txt', 'r').read().find('AGACATGCCT')
#TP53: 323910825

def motif_search(seq, subseq): 
	"""Search for a DNA motif in sequence. 
	use ambiguous values (like N = A or T or C or G, R = A or G etc.) 
	searches only on forward strand 
	""" 
	pattern = 'CCGCGNGGNGGCAG' 
	for motif in subseq: 
		value = IUPACData.ambiguous_dna_values[motif] 
		if len(value) == 1: 
			pattern += value 
		else: 
			pattern += '[%s]' % value 

	pos = -1 
	result = [pattern] 
	l = len(seq) 
	while True: 
		pos += 1 
		s = seq[pos:] 
		m = re.search(pattern, s) 
		if not m: 
			break 
		pos += int(m.start(0)) 
		result.append(pos) 
	return result

print(len(result))


#For giving updates to the user on the progress of the program
sys.stdeer.write()

#Output gene names, that contain binding motifs


#Import matrix as a numpy array???
#SCORING: The bigger the number the higher the preference
#>TP53_frequency_matrix
#A=[7544,10514,45,11931,1710,8,244,1228,1145,3851,3584,7104,0,12925,1825,327,879,1472]
#C=[1037,116,19689,204,374,2,12014,17286,11875,1474,756,26,19350,250,118,33,3338,9183]
#G=[6563,2433,13,653,137,20393,109,417,2151,9954,17846,8821,4,157,89,20549,116,825]
#T=[2268,735,837,739,18893,0,7106,2642,3774,2118,1546,190,1,1155,17864,68,17859,9120]

#>CTCF_frequency_matrix
#A=[87,167,281,56,8,744,40,107,851,5,333,54,12,56,104,372,82,117,402]
#C=[291,145,49,800,903,13,528,433,11,0,3,12,0,8,733,13,482,322,181]
#G=[76,414,449,21,0,65,334,48,32,903,566,504,890,775,5,507,307,73,266]
#T=[459,187,134,36,2,91,11,324,18,3,9,341,8,71,67,17,37,396,59]

#Assign scores 1-4, based on numbers in the frequency matrix
