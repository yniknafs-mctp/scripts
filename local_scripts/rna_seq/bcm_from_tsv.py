'''
Created on Jul 29, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np


'''
take tsv file and convert it to a big count matrix
'''

from ssea.lib.countdata_ysn import BigCountMatrix

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv_file")
    parser.add_argument("output_dir")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Creating a big count matrix from tsv file')
    
    bm = BigCountMatrix.from_tsv(args.tsv_file, args.output_dir)
    bm.close()
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
