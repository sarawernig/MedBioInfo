#! /usr/bin/env python3

from Bio import SeqIO	
from Bio.Alphabet import IUPAC		
from Bio.Seq import Seq
from Bio import motifs			
from Bio import SeqUtils

fn = open("input/frequency_matrixes_vertebrates_nr_JASPAR.txt","r")

with open("input/degenerateConsensusSequence.txt","a") as f:
	for m in motifs.parse(fn,"jaspar"):
		f.write(">" + str(m.base_id) + "\n")
		f.write(str(m.degenerate_consensus) + "\n")
		
fn.close()
