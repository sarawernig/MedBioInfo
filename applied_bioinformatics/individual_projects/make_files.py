#! /usr/bin/env python3

from Bio import SeqIO	#Import module for window sliding


#Change this to a string not an output file!
with open("windslid_10bp_out.txt","w") as f:
	for seq_record in SeqIO.parse("gencode.v26.lncRNA_transcripts.fa", "fasta"):
		for i in range(len(seq_record.seq) - 9) :		#Length of sequence minus 9
			f.write(">" + str(seq_record.id) + "\n")		#Write the name of the sequence
			f.write(str(seq_record.seq[i:i+10]) + "\n") 	#Window sliding for 10bp motifs, moving by 1bp	

with open("windslid_14bp_out.txt","w") as f:
	for seq_record in SeqIO.parse("gencode.v26.lncRNA_transcripts.fa", "fasta"):
		for i in range(len(seq_record.seq) - 13) :		#Length of sequence minus 13
			f.write(">" + str(seq_record.id) + "\n")		#Write the name of the sequence
			f.write(str(seq_record.seq[i:i+14]) + "\n") 	#Window sliding for 10bp motifs, moving by 1bp	


