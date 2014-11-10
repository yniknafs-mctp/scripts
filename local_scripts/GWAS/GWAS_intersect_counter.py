'''
Created on Jan 25, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np





    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("intersect_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    #count number of SNPs caught
    snps = set()
    for line in open(args.intersect_file):
        line = line.strip().split('\t')
        rsID = line[3]
        snps.add(rsID)
    for snp in snps: 
        print snp
    
    logging.info('Running main script')
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
