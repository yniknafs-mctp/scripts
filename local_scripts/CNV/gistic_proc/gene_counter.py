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
import collections


def function():
    return stuff

def main():
    # parse command line
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser()
    parser.add_argument("trans_file")
    parser.add_argument("peak_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    peak_file = open(args.peak_file)
    peak_headers = peak_file.next().strip().split('\t')
    q_dict = collections.defaultdict(lambda: 0)
    
    for line in peak_file: 
        line = line.split('\t')
        q_val = line[peak_headers.index('Residual q values after removing segments shared with higher peaks')]
        
        peakID = line[peak_headers.index('Unique Name')]
        q_dict[peakID] = q_val
    
    trans_file = open(args.trans_file)
    headers = trans_file.next().strip().split('\t')
    genes = []

    
    for line in trans_file: 
        line = line.split('\t')
        gene = line[headers.index('gene_id')]
        peakID = line[headers.index('peakID')]
        peak_type = peakID.split('Peak')[0]
        category = line[headers.index('category')]
        qval = q_dict[peakID]
        
        if (peak_type == "Amplification") and  (gene not in genes) and (category == 'intergenic'): 
            genes.append(gene)
    print len(genes)
        
    
    
    
    
    
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
