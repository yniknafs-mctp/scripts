from __future__ import division
import sys
import os
import copy
import string
import math
import time
import sets


def slice(list_id, col_id, col_data, filename_input):
	results = {}
	set_id = set(list_id)
	#read input file and place into list
	file_input = open(filename_input, 'r')
	for line in file_input:
		if line.strip():
			line_tabs = line.strip().split("\t")
			id = line_tabs[col_id].strip()
			if id in set_id:
				data = line_tabs[col_data].strip()
				results[id] = data
	file_input.close()
	
	#output data in same order as id list
	to_print = ""
	for id in list_id:
		if id not in results:
			to_print += "NA"
		else:
			to_print += results[id]
		to_print += "\t"
	
	return to_print



list_id = ["ENST00000395795","ENST00000374429","ENST00000395793","ENST00000374426","ENST00000343575","ENST00000395794"]
col_id = 0
col_data = 9

output_filename = "tcga_vs_breast_cxcl12.txt"

#list of all the samples that you want to run
sample_filename = "tcga_vs_breast.txt"

#paths to samples
path = "/tcga_vs/rnaseq_results/version_001_2013-01-14/"
#directory and name of bamfile in the sample directory
filesuffix = "/cufflinks_known/isoforms.fpkm_tracking"


#run code
sample_file = open(sample_filename, 'r')
output_file = open(output_filename, 'w')

#print header
header_str = "sample\t"
for id_str in list_id:
	header_str += id_str + "\t"
print >> output_file, header_str

i = 1
for sample_name in sample_file:
	sample_name = sample_name.strip()
	if sample_name:
		print >> sys.stderr, sample_name, "\tsample# ", i
		i += 1
		filename = path + sample_name + filesuffix
		result_str = slice(list_id, col_id, col_data, filename)
		print >> output_file, sample_name, "\t", result_str
		
sample_file.close()
output_file.close()





