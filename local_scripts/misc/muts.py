'''
Created on Jan 22, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections
import glob

'''
searches one maf file for mutations in a genes listed
'''

genes = [
          'HSD3B1', 
          'HSD3B2',
          'HSD17B1',
          'HSD17B2',
          'HSD17B3',
          'HSD17B4',
          'HSD17B5',
          'HSD17B6',
          'HSD17B7',
          'HSD17B8',
          'HSD17B9',
          'HSD17B10',
          'HSD17B11',
          'HSD17B12',
          'HSD17B13', 
          'HSD17B14',
          'UGT2B15', 
          'UGT2B17',
          'CYP17A1',
          'AKR1C1', 
          'AKR1C2',
          'AKR1C3',
          'AKR1C4',
          'CYP11A1',
          'CYP11B1',
          'CYP11B2',
          'CYP21A2',
          'CYP19A1',
          'SULT2A1',
          'AMFR',
          ]

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("mut_dir")
    parser.add_argument("gene")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    files = glob.glob(args.mut_dir + '/*')
    for fileh in files: 
        logging.debug('reading file: %s' % fileh)
        file = open(fileh)
        for line in file: 
            if line.startswith('#'):
                continue
            elif line.startswith("Hugo"): 
                header = line.strip().split('\t')
            line_s = line.strip().split("\t")
            if line_s[0] == args.gene:
                print fileh
                print '\t'.join(header)
                print line
    
if __name__ == '__main__': 
    sys.exit(main())
