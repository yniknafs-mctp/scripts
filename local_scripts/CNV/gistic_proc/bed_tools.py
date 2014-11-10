'''
Created on October 22, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv
import subprocess
from collections import defaultdict

logging.basicConfig(level=logging.DEBUG, 
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

peakID_handle = 3
transcriptID_handle = 7

def bed_handler(headers, lesion_file, compendia):
    logging.info("Generating bed file for intersection")
    
    bed_tmp = open('tmp.bed', 'w')
    
    genes_in_peak = defaultdict(lambda: [])
    trans_in_peak = defaultdict(lambda: [])
    for line in lesion_file: 
        wide_peak = line[headers.index("Wide Peak Limits")]
        wide_peak_list = wide_peak.replace('(',':').replace(')', '').replace('-',':').split(':')
        peak_loc = [wide_peak_list[0], wide_peak_list[1], wide_peak_list[2]]
        peak_name = [line[headers.index("Unique Name")]]
        peak = peak_loc + peak_name
        bed_tmp.write('\t'.join(map(str, peak)) + '\n')
    bed_tmp.close()
    
    logging.info('performing BEDTOOLS intersection')
    p1 = subprocess.check_output(["bedtools", 'intersect', '-a', 'tmp.bed', '-b', compendia, '-loj'])
    os.remove('tmp.bed')
    for line in p1.splitlines():
        line = line.split('\t')
        peakID = line[peakID_handle]
        if line[4] == '.':
            continue
        transcript_info= line[transcriptID_handle].split('|')
        transcriptID = transcript_info[1]
        geneID = transcript_info[0]
        if transcriptID not in trans_in_peak[peakID]:
            trans_in_peak[peakID].append(transcriptID)
        if geneID not in genes_in_peak[peakID]:
            genes_in_peak[peakID].append(geneID)
    return trans_in_peak, genes_in_peak
    
        
        