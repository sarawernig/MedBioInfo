#!/bin/bash

for filename in /home/sara/applied_bioinformatics/assignment1/yeast_genes/*.fasta

do

/home/sara/applied_bioinformatics/assignment1/muscle3.8.31_i86linux64 -in $filename -out "$filename"_muscle.fasta -log "$filename"_muscle.log &

done
