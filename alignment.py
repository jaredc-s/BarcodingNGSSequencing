# In his name
# Written by AREZOO MOVAGHAR
# Contact info: amovaghar@wisc.edu
#               muvaghar@gmail.com

# The propose of this code is aligning two sequences
# This program find all the possible alignments

import csv
import string
import sys
import time
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62

gap_open = -10

gap_extend = -0.5


def alignSeq(seq1,seq2,ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq, indels):
	alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)

	top_aln = alns[0]
	[aln1, aln2, score, begin, end] = top_aln
	ins_true=0
	ins_no=0
	del_true=0
	del_no=0
	fseq1 = []
	fseq2 = []
	res = []  # final result 
	n=0	
	shiftF=0
	if aln1.find('-')>0 or aln2.find('-')>0:
		indels[1]=indels[1]+freq
	
		
	for i in range(0, len(aln1)):
	
		if (aln1[i] == aln2[i]):
			n+=1
			mat[n]+=freq
			if ins_true==1:            
				in_freq[ins_no]+=freq  
				if (ins_no%3)==0:
					sh_in[n]+=freq; shiftF=1
				del_true=0
				ins_true=0
				ins_no=0
				del_no=0
			if del_true==1:
				del_freq[del_no]+=freq	
				if (del_no%3)==0:
						sh_del[n]+=freq; shiftF=1
				del_true=0
				ins_true=0
				ins_no=0
				del_no=0
		if aln1[i]=="-":	
				ins_true=1
				ins_no+=1
				ins[n]+=freq
				del_no=0
						
		if aln2[i]=="-":
				n+=1
				del_true=1
				del_no+=1
				dell[n]+=freq
				ins_no=0
		if aln2[i]!=aln1[i] and aln1[i]!="-" and aln2[i]!="-":		
				n+=1				
				mis[n]+=freq
				if ins_true==1:            
						in_freq[ins_no]+=freq   
						if (ins_no%3)==0:
							sh_in[n]+=freq; shiftF=1        	
				if del_true==1:           
						del_freq[del_no]+=freq
						if (del_no%3)==0:
							sh_del[n]+=freq; shiftF=1        	
		            
				del_true=0
				ins_true=0
				ins_no=0
				del_no=0	
		i+=1
	return shiftF

############################################################


def main(firstInput,wtInput,freqInput,res,freq_out):
	
	f = open(wtInput)
	seq1 = f.readline().upper()
	f.close()

	f2 = open(firstInput)
	f4 = open(freqInput)
#	print len(seq1)
	i=0
	lenSeq=len(seq1)+1
	ins=[0]*lenSeq
	dell=[0]*lenSeq
	mis=[0]*lenSeq
	mat=[0]*lenSeq
	sh_in=[0]*lenSeq
	sh_del=[0]*lenSeq
	in_freq=[0]*abs(lenSeq*2/3)
	del_freq=[0]*abs(lenSeq*2/3)
	indels=[0]*5
	for line in f2:
		freq=int(f4.readline())
		i+=1
		if line=="": break
		indels[0]=indels[0]+freq
		#print i
		seq2= line		
		#Call alignment function 
		shift=0
		shift=alignSeq(seq1.strip(),seq2.strip(),ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq, indels)
		if shift==1:
			indels[2]=indels[2]+freq
			#print indels[2]
		#print "****************"
	f2.close()
	f3=res
	f5=freq_out
	f6 = "indels.csv"
#print freq_out
	c1 = csv.writer(open(f3, "wt"))
	c2 = csv.writer(open(f5, "wt"))
	c3 = csv.writer(open(f6, "wt"))
	
	c1.writerow(' '+(seq1))
	c1.writerow(['Insert']+ins[1:])
	c1.writerow(['Delete']+dell[1:])
	c1.writerow(['Mismatch']+mis[1:])
	c1.writerow(['Match']+mat[1:])
	c1.writerow(['Insert_Inframe']+sh_in[1:])
	c1.writerow(['Delete_Inframe']+sh_del[1:])
	c2.writerow(['length']+range(1,abs(lenSeq*2/3)))
	c2.writerow(['Insert']+in_freq[1:])
	c2.writerow(['Delete']+del_freq[1:])
	indels[3]=indels[0]-indels[1] #wt
	indels[4]=indels[1]-indels[2] #indels without frameshift

	c3.writerow(' '+ ['Total']+['Indels']+['FrameShift']+	['WT+	Indels-Frameshifts'])
	#print indels
	c3.writerow(['indels']+indels)
	#print "end"
	return 0
############################################################

main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
