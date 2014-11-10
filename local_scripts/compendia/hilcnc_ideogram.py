'''
Created on Apr 24, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("meta")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    clincs = []
    for line in meta_fh: 
        line = line.strip().split('\t')
        uce = line[meta_header.index('uce')]
        chrom = line[meta_header.index('chrom')]
        start = line[meta_header.index('start')]
        end = line[meta_header.index('end')]
        strand = line[meta_header.index('strand')]
        ideo = line[meta_header.index('ideo_cat')]
        if ideo == 'cscat':
            lineo = ['mapping', '6', 'red', '.', 'box', chrom, start, end, strand]
            print '\t'.join(lineo)
        if ideo == 'ls':
            lineo = ['mapping', '6', 'blue', '.', 'box', chrom, start, end, strand]
            print '\t'.join(lineo)
#         
        
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

