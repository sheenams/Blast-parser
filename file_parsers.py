# ~/file_parsers.py
#!/usr/bin/env python
import os, sys, csv, re, gzip
from os import path
from collections import Counter

from Parser import *

try: 
	from Bio import SeqIO
	from Bio.Blast import NCBIXML

except ImportError:
    sys.stderr.write("BioPython not installed correctly (see http://biopython.org)\n")
    sys.exit(-1)
    

#Set parameters for filtering
IDENTITY = results.IDENTITY

LENGTH= results.LENGTH

E_VALUE= results.E_VALUE

BIT_SCORE= results.BIT_SCORE

print >> sys.stderr, '''Parsing parameters set to: %s identity, %s length, %s e_value, %s bit_score. 
If these are not correct, stop the script and set with -i, -l, -e or -b respectively \n''' %(IDENTITY,LENGTH,E_VALUE, BIT_SCORE )


##Parses the database file (ends with fna) and returns the genome list
##Exists program if database file not found
def parse_db_file(filenames):
	genome_list=[]
	all_genomes=[]
	fna_file_object=None
	for filename in filenames:
		if os.path.isfile(filename) and filename.endswith('.fna'):
			##open the file and store in filename_object
			fna_file_object= open(filename, 'r')
			for seq_record in SeqIO.parse(fna_file_object, 'fasta'):
				genome=seq_record.id.split()
				genome= genome[0].split('_')
				genome=genome[0]
				genome=str(genome)
				all_genomes.append(genome)
				##make list of all unique genomes (global) 
				genome_list=set(all_genomes)
				#genome_list=sorted(genome_list)
			print >> sys.stderr, 'DATABASE FILE FOUND:', filename
	if not fna_file_object:
		print >> sys.stderr, 'No database file found.\n'
		sys.exit(1) 
			
	return genome_list


##Prints genome header for output file
##Provides genome list for blast output parsing 
def print_genome_header(genome_list, flag):
	if not flag:
		for entry in genome_list:
			print '\t', entry, 
		print '\t'
		flag = True
	return flag
	
	
	
##Uses genome list from database file to print output
##Parameters provided at runtime or uses all defaults of 0
##Parse xml file 	
def parse_xml_file(filename, genome_list):
	##open files
	if filename.endswith(".gz"):
		result_handle=gzip.open(filename, 'r')
		print >> sys.stderr, 'GZ FILE FOUND:', filename
			
	elif filename.endswith(".out"):
		result_handle= open(filename, 'r')
		print >> sys.stderr, 'XML FILE FOUND:', filename
	
	##actual parsing: 	
	BLAST_RECORDS = NCBIXML.parse(result_handle)	
	##Create Q dictionary 
	Q=dict()
	##Go through each BLAST_RECORD, save each QUERY/GENOME that fits the parameters
	for ONE_RECORD in BLAST_RECORDS:
		HEADER = ONE_RECORD.query
		HEADER = HEADER.split()
		QUERY = HEADER[0]
		QUERY = str(QUERY)
		##fill Q dictionary with QUERY dictionary (list of query numbers)
		if not QUERY in Q:
			Q[QUERY]={}	
		for alignment in ONE_RECORD.alignments:
			SUBJECT_ID = alignment.title
			SUBJECT_ID = SUBJECT_ID.split()
			GENOME= SUBJECT_ID[1].split('_')
			GENOME= GENOME[0]
			GENOME = str(GENOME)
			for hsp in alignment.hsps:
				PERCENT_IDENTITY= round((100.0*(hsp.identities)/ (hsp.align_length)),2)
				if (PERCENT_IDENTITY >= IDENTITY and alignment.length >= LENGTH and hsp.expect >= E_VALUE and hsp.bits >= BIT_SCORE):
					##fill QUERY dictionary with GENOME keys
					if not GENOME in Q[QUERY].keys():
						Q[QUERY][GENOME]=1
																
	##use Counter module to count the GENOMES present in the QUERY dictionary
	genome_count=Counter()
	for QUERY in Q.keys():
		for GENOME in Q[QUERY].keys():
			genome_count[GENOME] +=1

	##Find total amount of queries in file, includes queries with no hits 
	query_length=len(Q)
	print 'Number of queries: ', query_length
	
	##Print results that met parameters where:
	##Any amount of genome hits per query = 1
	##Output ratio of genome hits divided by total queries
	line = str(filename)
	for GENOME in genome_list:
		out = "%.4f" % genome_count[GENOME]
		out = float(out)/query_length
		out = "%.4f" % out
		line = line + '\t' + str(out)
	print line
	print >> sys.stderr, 'XML file parsed: ',filename

			
#Parse blast table file 	

##Uses genome list from database file to print output
##Parameters provided at runtime or uses all defaults of 0
##Parse xml file
def parse_blast_tabfile(filename, genome_list):
	##open files
	with open(filename, 'r') as f:
		print >> sys.stderr, 'TABLE FILE FOUND: ', filename
		Q = {}
		reader = csv.reader(f, delimiter="\t")
		for row in reader:
			##matches lines that don't start with a number
			comment_line=re.match('[^0-9]', row[0]) 
			if comment_line:
				##matches the # Query lines
				##Ensures all Queries are accounted for, even if it had no hits
				query_line = re.match('#\sQ.*', row[0])
				if query_line:
					QUERY= row[0].split(' ')
					QUERY=QUERY[2]
					if not Q.has_key(QUERY):
						Q[QUERY]={}
	
			if not comment_line:
				QUERY=str(row[0])
				if not QUERY in Q.keys():
					Q[QUERY]={}
				GENOME=row[1].split('_')
				GENOME=GENOME[0]
				row_identity = float(row[2])
				row_length= str(row[3])
				row_e_value=float(row[10])
				row_bit_score=float(row[11])
				if (row_identity >= IDENTITY and row_length >=LENGTH and row_e_value >=E_VALUE and row_bit_score >=BIT_SCORE):
					if not GENOME in Q[QUERY].keys():
						Q[QUERY][GENOME]=1
															
	###use Counter module to count the GENOMES present in the QUERY dictionary
	genome_count=Counter()
	for QUERY in Q.keys():
		for GENOME in Q[QUERY].keys():
			genome_count[GENOME] +=1

	###Find total amount of queries in file
	query_length=len(Q)
	print 'Number of queries: ', query_length

	#Print results that met parameters where:
	#Any genome hit per query = 1
	#Output ratio of genome hits divided by total queries, even if query had no hits
	line = str(filename)
	for GENOME in genome_list:
		out = "%.4f" % genome_count[GENOME]
		out = float(out)/query_length
		out = "%.4f" % out
		line = line + '\t' + str(out)
	print line
	print >> sys.stderr,'TABLE file parsed:', filename
	
		
			
			
			
