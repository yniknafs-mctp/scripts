'''
Created on Jun 11, 2014
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
    parser.add_argument("metadata")
    parser.add_argument("seq")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Saving sequences')
    seq_dict = {}
    seq_fh = open(args.seq)
    for line in seq_fh:
        line = line.strip().split('\t')
        t_id = line[0].split('|')[1]
        seq = line[2].upper()
        seq_split = seq.split('|')
        seq_out = '<->'.join(seq_split)
        ''
        seq_dict[t_id] = seq_out
    
    meta_fh = open(args.metadata)
    meta_header = meta_fh.next().strip().split('\t')
    meta_header.append('seq')
    print '\t'.join(meta_header)
    for line in meta_fh:
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        seq = seq_dict[t_id]
        line.append(seq)
        print '\t'.join(line)
        
    
    
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
