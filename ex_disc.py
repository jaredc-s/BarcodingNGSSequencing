import collections
import sys
import csv

#CDH20 ,atccagacagtgagtgcggt,	cctcaggagccaagctgtag, ctacagcttggctcctgagg		
#CDH5 ,gactccttccagcttcacca,	cacggacgcattgaacaacc, ggttgttcaatgcgtccgtg
#FMR1 ,TGCCACCAAATTCCCTTCCT,	AGGTTGAGAAAAATGGGTGCT, AGCACCCATTTTTCTCAACCT
#FMR1 ,TGCCACCAAATTCCCTTCCT,	AGGTTGAGAAAAATGGGTGCT, AGCACCCATTTTTC reduced
#FOXP2 ,agcatctgctcagccttcag,	ttgaggcagcgattggacag, ctgtccaatcgctgcctcaa
#GAPDH ,acccagaagactgtggatg,	gtagaggcagggatgatgt,acatcatccctgcctctac		
#NANOG ,aaggcctcagcacctaccta,	GAAGGTTCCCAGTCGGGTTC,GAACCCGACTGGGAACCTTC		
#OCT4 ,gtggtccgagtgtggtt,	gaaaggagacccagcag, ctgctgggtctcctttc
#Pax6 ,aagcaaaatagcccagtataag,tatgttatcgttggtacagacc,ggtctgtaccaacgataacata
#POU4F1 ,cctcacccgcatatgcacag,CCCGGACGGCATGTTCA, TGAACATGCCGTCCGGG
	

def main(SeqInput,SeqOutput,freqOutput):
	sub=("CATAGCCGTATAG").upper()
	sub2=("CACGTCTGAACTC").upper()
	
	f = open(SeqInput)
	
	c1 = open(SeqOutput, "wt")
	c2 = open(freqOutput, "wt")
	#c3 = open(wt, "wt")
	words = collections.Counter(f)
	l= words.most_common(len(words))
#	strt_ind = l[0][0].upper().find(sub)
#	end_ind= l[0][0].upper().find(sub2)
#	
#	c3.write((l[0][0])[strt_ind:end_ind+len(sub2)])
	i=0	
	while i<len(l) and (l[i][1])>=100:
		strt_ind = l[i][0].upper().find(sub)
		end_ind= l[i][0].upper().find(sub2)
		if strt_ind!=-1 and end_ind!=-1:				
			c1.write((l[i][0])[strt_ind:end_ind+len(sub2)])
			c1.write('\n')			
			c2.write(str(l[i][1]))
			c2.write('\n')
		i+=1

	return 0
main(sys.argv[1], sys.argv[2], sys.argv[3])
