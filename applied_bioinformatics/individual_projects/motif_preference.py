#! /usr/bin/env python3

from Bio import SeqIO   #For window sliding


with open("windslid_out.txt","w") as f:
        for seq_record in SeqIO.parse("gencode.v26.lncRNA_transcripts.fa", "fasta"):
            for i in range(len(seq_record.seq) - 9) :		#Length of sequence minus 9
               f.write(str(seq_record.id) + "\n")		#Write the name of the sequence
               f.write(str(seq_record.seq[i:i+10]) + "\n") 	#Window sliding for 5bp motifs, moving by 1bp	

#Make a dictionary with TFs and their binding motifs
TP53=AGACATGCCT		#TP53 10bp binding motif

#The bigger the number the higher the preference
#>TP53_frequency_matrix
#A=[7544,10514,45,11931,1710,8,244,1228,1145,3851,3584,7104,0,12925,1825,327,879,1472]
#C=[1037,116,19689,204,374,2,12014,17286,11875,1474,756,26,19350,250,118,33,3338,9183]
#G=[6563,2433,13,653,137,20393,109,417,2151,9954,17846,8821,4,157,89,20549,116,825]
#T=[2268,735,837,739,18893,0,7106,2642,3774,2118,1546,190,1,1155,17864,68,17859,9120]


#Assign scores 1-4, based on numbers in the frequency matrix
