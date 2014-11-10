import argparse
import logging
import os
import sys
import csv
import time
import collections
import bisect
import operator
import subprocess

'''
Takes the tsv for SNPS (after collapsing by SNPid) and converts to BED

input: unique SNP tsv (gwas_rsID.tsv from GWAS_collapse_by_rsID.py)
'''

def main():

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("GWAS_file")
    args = parser.parse_args()
    logging.info("Starting")
    
    GWAS_BED = open("gwas.bed", 'w')
    with open(args.GWAS_file, 'r') as f:
        headers = f.next().strip().split('\t')
        for row in f:
            row = row.strip().split('\t')
            if row != ['']:
                ID = row[headers.index('Snp_id_current')]
                chrom = 'chr'+row[headers.index('Chr_id')]
                start = row[headers.index('Chr_pos')]                
                end = int(start) + 1
                line = [chrom, start, end, ID]
                print >> GWAS_BED, '\t'.join(map(str, line))
    GWAS_BED.close()
            
    
            
    return 0



if __name__ == '__main__':

    sys.exit(main())