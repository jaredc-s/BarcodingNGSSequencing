# Written by AREZOO MOVAGHAR
# Contact info: amovaghar@wisc.edu
#               muvaghar@gmail.com
# This script reads the fastq file, copies the main sequence to the new file.
# python readseq.py input.fastq output.txt
 
import collections
import sys
import csv

def main(SeqInput,SeqOutput):
	f = open(SeqInput)
	f2 = open(SeqOutput,'w+')
	sub=("CCTCACCCGCATATGCACAG").upper()
	sub2=("TGAACATGCCGTCCGGG").upper()
	while 1:
		s=f.readline()
		if not s:
			break
		l=f.readline()
		strt_ind = l.upper().find(sub)
		end_ind= l.upper().find(sub2)		

		if strt_ind!=-1 and end_ind!=-1:				
			f2.write(l[strt_ind:end_ind+len(sub2)])
			f2.write('\n')			

		f.readline()
		f.readline()
	f.close()
	f2.close()
	return 0
	
main(sys.argv[1], sys.argv[2])
	