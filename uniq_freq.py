import collections
import sys
import csv



# python uniq_freq.py Input_Seq.txt Output_listOfUniqueSeq.txt Output_Frequency.txt	

def main(SeqInput,SeqOutput,freqOutput):

	
	f = open(SeqInput)
	
	c1 = open(SeqOutput, "wt")
	c2 = open(freqOutput, "wt")
	words = collections.Counter(f)
	l= words.most_common(len(words))

	i=0	
	while i<len(l) and (l[i][1])>=100:				
			c1.write((l[i][0]))
			c2.write(str(l[i][1]))
			c2.write('\n')
			i+=1

	return 0
main(sys.argv[1], sys.argv[2], sys.argv[3])
