#! /usr/bin/env python3
import string
import random

length = int(input("Length:"))	
def random_name(size=10, chars=string.ascii_lowercase):
	return '>' + ''.join(random.choice(chars) for _ in range(size))

random_name()

def random_fasta(size=length, chars='GATC'):
	#''.join joins an empty string '' with randomly generated characters
	return ''.join(random.choice(chars) for _ in range(size))

random_fasta()

print (random_name(), '\n', random_fasta())



	



