'''
Created on Apr 15, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections


'''
read gwas metadata file and add info for the gwas snp 
'''


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("meta")
    parser.add_argument("gwas")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
     
    #read gwas file and store info regarding gwas snps
    gwas_fh = open(args.gwas)
    gwas_header = gwas_fh.next().strip().split('\t')
    gwas_canc_dict = collections.defaultdict(lambda: '0')
    for line in gwas_fh:
        line = line.strip().split('\t')
        snp = line[gwas_header.index('SNPs')]
        tissue = line[gwas_header.index('Tissue')]
        gwas_canc_dict[snp] = tissue
        
    
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_header.append('gwas_cancer')
    print '\t'.join(meta_header)
    for line in meta_fh:
        line = line.strip().split('\t')
        snp = line[meta_header.index('gwas_snp_id')]
        gwas_can = gwas_canc_dict[snp]
        line.append(gwas_can)
        print '\t'.join(line)
        
        
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
