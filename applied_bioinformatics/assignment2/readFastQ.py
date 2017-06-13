#! /usr/bin/env python3

import sys
import re

#If user provides a filename in the command line when running the program, than use that file, else ask the user to provide the filename.

##filename = 'example_fastQ_file.fastq'

if len(sys.argv) == 1:
	filename = input("Enter Filename: ")
else:
	filename = sys.argv[1]

#Count number of occurances of @ symbol
LineNumber = 0

InFile = open(filename, 'r')

for Line in InFile:
	Line = Line.strip('\n')
	LineNumber = LineNumber + 1

FastaCount = "%d" % (LineNumber / 4)
print('The file contains', FastaCount, 'FastQ sequences.' )

#Print the lines containing '@', without '@'

searchString = "(^[@])(\w+)"
Result = re.search(searchString, InFile)
Result.group(2)

print(Result.group(2))	

InFile.close()

##You are not allowed to read the file more than once!!!!
