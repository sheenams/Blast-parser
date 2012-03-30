# ~/blast_parse_final.py
import sys
#import csv
#import pickle

FILE_INPUT = open("blast_test.xml")
#FILE_OUTPUT=open('output_test.csv','w')

#Get parse parameters from user
print 'Please enter the following parameters. Enter 0 if not applicable'

IDENTITY = int(raw_input("Minimum identity: "))
LENGTH= int(raw_input("Minimum length: "))
E_VALUE= int(raw_input("Minimum e-value: "))
BIT_SCORE= int(raw_input("Minimum bit score: "))
#UNIQUE= raw_input("> ")

#Import XML parser from BioPython
from Bio.Blast import NCBIXML
sys.stdout=open('out_data_check2.txt', 'w')
#Iterate over entire file and return Blast Record for each query
BLAST_RECORDS = NCBIXML.parse(FILE_INPUT)

#Read over each record 
#ONE_RECORD = BLAST_RECORDS.next()
Q=dict()
genome_list=[]
for ONE_RECORD in BLAST_RECORDS:
	HEADER = ONE_RECORD.query
	HEADER = HEADER.split()
	QUERY = HEADER[0]
	QUERY = str(QUERY)
	for alignment in ONE_RECORD.alignments:
		SUBJECT_ID = alignment.title
		SUBJECT_ID = SUBJECT_ID.split()
		GENOME= SUBJECT_ID[1].split('_')
		GENOME= GENOME[0]
		GENOME = str(GENOME)
		for hsp in alignment.hsps:
			PERCENT_IDENTITY= round((100.0*(hsp.identities)/ (hsp.align_length)),2)
			if PERCENT_IDENTITY >= IDENTITY:
				if alignment.length >= LENGTH:
					if hsp.expect >= E_VALUE:
						if hsp.bits >= BIT_SCORE:
							#print '\n**Next Record** \n'
							#print 'Subject ID (Genome):', GENOME
							##fill Q dictionary with Query dictionary (list of query numbers)
							if not QUERY in Q:
								Q[QUERY]={}	
								#fill Query dictionary with Genome keys and count them
							if not GENOME in Q[QUERY].keys():
								Q[QUERY][GENOME]=1
							else:	
								Q[QUERY][GENOME]+=1
							genome_list.append(GENOME)

#remove duplicates from genome list
genome_list=set(genome_list)
print 'genome list:', genome_list

#For every query in the dictionary
for entry in genome_list:
	print 'Genome:', entry
	for QUERY in Q.keys():
	#check genome in list is in the query dictionary
	#if not, print '0'
		if not entry in Q[QUERY].keys():
			print 'Query:', QUERY, '0'
		else:
			print 'Query:', QUERY, '1'

#output=open('myfile.pkl','wb')
#pickle.dump(Q, output)
#output.close()


#writing = csv.DictWriter(FILE_OUTPUT, fieldnames='QUERY')
#fields = genome_list
#w=csv.DictWriter(fid, filed)
#w.writerow(dict(zip(fields, fieds)))
#for row in rows:
	#w.writerow(row)
