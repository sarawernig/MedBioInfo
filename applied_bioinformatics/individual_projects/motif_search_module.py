#! /usr/bin/env python3


import re		#Import module for regular expressions
from Bio import SeqIO	#Import module for window sliding
from Bio.Alphabet import IUPAC		#Import module for ambigious sequences
from Bio.Seq import Seq


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
