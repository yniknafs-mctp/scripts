'''
Created on October 28, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv
import collections

rna_id_file = '/home/yniknafs/swap/rna_seq/rnaseq_library_info.txt'


def id_converter(matrix_file):
#create DNA_id_file from clinicalMatrix TCGA data (i.e. match TCGA ptID to the CNV sampleID)
#used to convert CNV_DNA sampleID's to RNA sampleID's which are used in the expression file
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Creating ID conversion dictionary')
    pt_to_RNA = collections.defaultdict(lambda: [])
    fileh = open(rna_id_file)
    header = fileh.next().strip().split("\t")
    for line in fileh: 
        line = line.strip().split('\t')
        RNA_sample = line[header.index('library_id')]
        pt = line[header.index('tcga_legacy_sample_id')]
        if pt.startswith("TCGA"):
            pt_split = pt.split('-')
            pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3][:-1]])
            pt_to_RNA[pt_new] = RNA_sample
    
    DNA_to_RNA = collections.defaultdict(lambda: [])
    with open(matrix_file, 'r') as f:
        headers = f.next().strip().split('\t')
        print headers 
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            CNV_sampleID = row[0]
            ptID = row[1]
    #building dictionary to convert GISTIC id to library ID   
            DNA_to_RNA[CNV_sampleID] = pt_to_RNA[ptID][:]
    
    return DNA_to_RNA
