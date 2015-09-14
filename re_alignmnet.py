# Written by AREZOO MOVAGHAR
# Contact info: amovaghar@wisc.edu
#               muvaghar@gmail.com
# This script reads the list of unique sequences, wiletype sequence and frequency of sequences as the input, perfurms alignment and retuens two output files.
# The first output file lists frequency of  insertions,deletions,mismathes, matches, inframe insertions and inframe deletions based on the position.
# The second output file lists frequency of insertions and deletions based on the length
# python re_alignment Input_listOfUniqueSeq.txt Input_wildtype.txt Input_Frequency.txt Output_Alignments_position.txt Output_Alignment_length.txt

import csv
import string
import sys
import time
# Create score matrix and direction matrix (to handle trace back) 
def init(seq1,seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    
	#linear gap penalty and scores 
    d=3 # gap penalty
    matchScore=2
    misScore=1 
    
    # Create score matrix 
    F = [[0 for col in range(len1+1)] for row in range(len2+1)]

    # Create dir matrix ( will be used in back track part)
    dir = [[0 for col in range(len1+1)] for row in range(len2+1)]
    for row in range(0,len2+1):
        for col in range(0,len1+1):
            dir[row][col] = []

	
    # Initialize the first row and column using the gap penalty value (d=1)
    for col in range(1,len1+1):
        F[0][col] = F[0][col-1] - d
        dir[0][col] = ['W']

    for row in range(1,len2+1):
        F[row][0] = F[row-1][0] - d
        dir[row][0] = ['N']

    # Calculate other values 
    # N means North (i-1) and valN means f-value from i-1
    # W means West (j-1) and valW means f-value from j-1 
    # mat and mis refer to i-1 j-1. We can use just NW  and the final result will be fine but considering them seprately will improve our results and gaps will be okay for seq1 too ;)
    
    for row in range(1,len2+1):
        for col in range(1,len1+1):
            if (seq1[col-1] == seq2[row-1]):
                valMat = F[row-1][col-1] + matchScore
                temp = 'Mat'
            else:
                valMat = F[row-1][col-1]+misScore
                temp = 'Mis'
            valN = F[row-1][col] - d
            valW = F[row][col-1] - d
            value = max(valMat,valN,valW)
            F[row][col] = value
            if (value == valMat):
                dir[row][col].append(temp)
            if (value == valN):
                dir[row][col].append('N')
            if (value == valW):
                dir[row][col].append('W')
    return [F, dir]
############################################################

# finding all alignments using direction matrix recursively
# keep the result in path 
def alignFunc(dir, row, col, t):
    waiting=time.time() -t
    #print waiting
    if waiting> 0.05:
		return [[]]

    if (row == 0) and (col == 0): # return 
        return [[]]
    paths = []

    if ('Mat') in dir[row][col]:  # recursion to i-1, j-1 
        for rec in alignFunc(dir, row-1, col-1,t):
            paths.append(rec + ['Mat'])
            
    if ('Mis') in dir[row][col]: # recursion to i-1, j-1 
        for rec in alignFunc(dir, row-1, col-1,t):
            paths.append(rec + ['Mis'])

    if ('W') in dir[row][col]:   # recursion to i-1, j
        for rec in alignFunc(dir, row, col-1,t):
            paths.append(rec + ['W'])
            
    if ('N') in dir[row][col]:  # recursion to i, j-1
        for rec in alignFunc(dir, row-1, col,t):
            paths.append(rec + ['N'])
            
    return paths 

############################################################
# The final result 
def traceback(alignPath, seq1, seq2,ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq):

    ind1 = 0  # pointer to current position in seq1
    ind2 = 0  # pointer to current position in seq2
    ins_true=0
    ins_no=0
    del_true=0
    del_no=0
    fseq1 = []
    fseq2 = []
    res = []  # final result 
    n=0
    # Build the alignment sequences based on the alignment path     
    # I print aligend seq1, aligned seq2 and the final alignment because sometimes - should apear in the first seq! 
    for step in alignPath: 
        
        if (step == 'N'):  # j-1  
            fseq1.append('-')
            ins[n]+=freq
            fseq2.append(seq2[ind2])
#            res.append('-')
            ins_no+=1
            ins_true=1;            
            ind2 = ind2 + 1
         	
        elif (step == 'W'): # i-1
            n+=1
            #print "del"        
            fseq1.append(seq1[ind1])
            fseq2.append('-')
            dell[n]+=freq            
        #    res.append('-')
            del_no+=1
            del_true=1
                      	

            ind1 = ind1 + 1
            
        elif (step == 'Mat'):
            #print "Mat"  
            n+=1      
            fseq1.append(seq1[ind1])            
            fseq2.append(seq2[ind2])
            mat[n]+=freq
       #     res.append(seq1[ind1])       		
            if ins_true==1:            
				in_freq[ins_no]+=freq  
				if (ins_no%3)==0:
					sh_in[n]+=freq          	
            if del_true==1:           
				del_freq[del_no]+=freq	
				if (del_no%3)==0:
						sh_del[n]+=freq             						
            
            del_true=0
            ins_true=0
            ins_no=0
            del_n=0
            ind1 = ind1 + 1
            ind2 = ind2 + 1
            mis[n]+=freq
        elif (step == 'Mis'):  
           # print "Mis"
            n+=1        
            fseq1.append(seq1[ind1])
            fseq2.append(seq2[ind2])
                     
            if ins_true==1:            
				in_freq[ins_no]+=freq  
				if (ins_no%3)==0:
					sh_in[n]+=freq          	
            if del_true==1:           
				del_freq[del_no]+=freq	
				if (del_no%3)==0:
						sh_del[n]+=freq             
            
            del_true=0
            ins_true=0
            ins_no=0
            del_n=0     
            ind1 = ind1 + 1
            ind2 = ind2 + 1    

        else: # defualt short story 
            fseq1.append('?')
            fseq2.append('!')
     #       res.append('O')
    if ins_true==1:
         	in_freq[ins_no]+=freq 
           	#print "end,",ins_no          	
    if del_true==1:           
			del_freq[del_no]+=freq
			#print "endd",del_no
    s1=string.join(fseq1, '')
    s2=string.join(fseq2, '')
    return [s1,s2]
############################################################

def alignSeq(seq1,seq2,ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq):
    (vals,dir) = init(seq1,seq2)
    a=""

    t=time.time()
    alignPath=alignFunc(dir, len(seq2), len(seq1),t)[0]
    print alignPath

    [a,b]=traceback(alignPath, seq1, seq2,ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq)  # Based on the path we show the results
    print a
    print b

        
    return [a,b]
############################################################

def main(uniqueSeqFile,wtSeqFile,freqInput,res,freq_out):
	f = open(wtSeqFile)
	wtSeq = f.readline()
	f.close()

	f2 = open(uniqueSeqFile)
	f4 = open(freqInput)

	i=0
	ins=[0]*len(wtSeq)
	dell=[0]*len(wtSeq)
	mis=[0]*len(wtSeq)
	mat=[0]*len(wtSeq)
	sh_in=[0]*len(wtSeq)
	sh_del=[0]*len(wtSeq)
	in_freq=[0]*21
	del_freq=[0]*21
	for line in f2:
		freq=int(f4.readline())
		i+=1
		if line=="": break
		#print i
		seq2= line
		#Call alignment function 
		[a,b]=alignSeq(wtSeq.strip(),seq2.strip(),ins,dell,mis,mat,sh_in,sh_del,in_freq,del_freq,freq)
		#print "****************"
	f2.close()
	f3=res
	f5=freq_out
	print freq_out
	c1 = csv.writer(open(f3, "wt"))
	c2 = csv.writer(open(f5, "wt"))
	c1.writerow(' '+(wtSeq))
	c1.writerow(['Insert']+ins[1:])
	c1.writerow(['Delete']+dell[1:])
	c1.writerow(['Mismatch']+mis[1:])
	c1.writerow(['Match']+mat[1:])
	c1.writerow(['Insert_Inframe']+sh_in[1:])
	c1.writerow(['Delete_Inframe']+sh_del[1:])
	c2.writerow(['length']+range(1,21))
	c2.writerow(['Insert']+in_freq[1:])
	c2.writerow(['Delete']+del_freq[1:])

	#print "end"
	return 0
############################################################

main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
