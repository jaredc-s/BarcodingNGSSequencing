# Written by AREZOO MOVAGHAR
# Contact info: amovaghar@wisc.edu
#               muvaghar@gmail.com
# This script reads the input sequence file, makes an index of the sequences and returns a list of unique sequences and frequency of each sequence in the input file. It also returns the most frequent sequence as the wildtype. 
# python dicts.py Input_Seq.txt Output_listOfUniqueSeq.txt Output_Frequency.txt Output_wildtype.txt

import collections
import sys
import csv

def main(SeqInput,SeqOutput,freqOutput,wt):
	f = open(SeqInput)
	
	c1 = open(SeqOutput, "wt")
	c2 = open(freqOutput, "wt")
	c3 = open(wt, "wt")
	words = collections.Counter(f)
	l= words.most_common(len(words))
	c3.write(l[0][0])

	for i in range (0, len(l)):
		print l[i]
		c1.write(l[i][0])
		c2.write(str(l[i][1]))
        c2.write('\n')
    

	return 0
main(sys.argv[1], sys.argv[2], sys.argv[3],sys.argv[4])
