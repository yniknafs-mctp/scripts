'''
Created on Apr 3, 2014
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
    parser.add_argument("tips")
    parser.add_argument("ssea")
    parser.add_argument("dg")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    tip_dict = collections.defaultdict(lambda: set())
    for line in open(args.tips):
        line = line.strip().split('\t')
        tip_dict[line[1]].add(line[0])
    
    
    
    ssea_fh = open(args.ssea)
    ssea_header = ssea_fh.next().strip().split('\t')
    print '\t'.join(ssea_header)
    for line in ssea_fh: 
        line = line.strip().split('\t')
        peak = line[ssea_header.index('ss_compname')]
        t_id = line[ssea_header.index('transcript_id')]
        if peak == 'cnv_cancer_vs_normal':
            print '\t'.join(line)
        if t_id in tip_dict[peak]:
            print '\t'.join(line)
        
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
