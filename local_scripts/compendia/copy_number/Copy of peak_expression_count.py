'''
Created on Apr 6, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections


LAST_PHENO_COL = 13

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("peak_smx")
    parser.add_argument("ca_type_smx")
    parser.add_argument("expr")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        
    logging.info('fuck')
    peak_fh = open(args.peak_smx)
    peak_fh.next()
    peak_fh.next()
    peak_dict = collections.defaultdict(lambda: 'mismatch')
    for line in peak_fh: 
        line = line.strip().split('\t')
        pt = line[0]
        val = line[1]
        peak_dict[pt] = val
    
    type_fh = open(args.ca_type_smx)
    type_fh.next()
    type_fh.next()
    type_dict = collections.defaultdict(lambda: 'mismatch')
    for line in type_fh: 
        line = line.strip().split('\t')
        pt = line[0]
        val = line[1]
        type_dict[pt] = val
    
    expr_fh = open(args.expr)
    expr_head = expr_fh.next().strip().split('\t')
    
    headero = [
               'transcript_id',
               'gene_id',
               'pt_id',
               'cnv_val',
               'type_val',
               'expr_val'
               ]
    print '\t'.join(headero)
    for line in expr_fh: 
        line = line.strip().split('\t')
        t_id = line[expr_head.index('tracking_id')]
        g_id = line[expr_head.index('gene_id')]
        for x in xrange(LAST_PHENO_COL, len(line)):
            expr_val = line[x]
            pt = expr_head[x]
            peak_val = peak_dict[pt]
            type_val = type_dict[pt]
            lineo = t_id, g_id, pt, peak_val, type_val, expr_val
            print '\t'.join(map(str, lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
