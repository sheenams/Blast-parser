# ~/blast_parser_final.py
#!/usr/bin/env python

##Import modules needed for code
import sys, os, gzip
import file_parsers

flag=False

filenames = os.listdir(os.curdir)


genome_lists=[]

genome_list= file_parsers.parse_db_file(filenames)

file_parsers.print_genome_header(genome_list, flag)



for filename in filenames:
	
	if os.path.isfile(filename) and filename.endswith(".fna"):
		continue
	
	elif os.path.isfile(filename) and (filename.endswith("xml.blast.out") or filename.endswith("xml.blast.out.gz")):
		file_parsers.parse_xml_file(filename, genome_list)  
	
	elif os.path.isfile(filename) and filename.endswith("blast.out"):
		file_parsers.parse_blast_tabfile(filename, genome_list) 

