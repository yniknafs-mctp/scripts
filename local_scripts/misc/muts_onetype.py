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
searches directories with mutation calls for a given gene and identifies all mutations for that gene
for TCGA
'''

fields = [
          'Hugo_Symbol',
          'file',
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
