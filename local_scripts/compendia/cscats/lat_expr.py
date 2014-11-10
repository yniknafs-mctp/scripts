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


LAST_PHENO_COL = 1

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("ca_vs_norm_smx")
    parser.add_argument("ca_type_smx")
    parser.add_argument("pancan_smx")
    parser.add_argument("expr")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        
    logging.info('fuck')
    canc_fh = open(args.ca_vs_norm_smx)
    canc_fh.next()
    canc_fh.next()
    canc_dict = collections.defaultdict(lambda: 'mismatch')
    for line in canc_fh: 
        line = line.strip().split('\t')
        pt = line[0]
        val = line[1]
        canc_dict[pt] = val
    
    type_fh = open(args.ca_type_smx)
    type_fh.next()
    type_fh.next()
    type_dict = collections.defaultdict(lambda: 'mismatch')
    for line in type_fh: 
        line = line.strip().split('\t')
        pt = line[0]
        val = line[1]
        type_dict[pt] = val
    
    pancan_fh = open(args.pancan_smx)
    pancan_fh.next()
    pancan_fh.next()
    pancan_dict = collections.defaultdict(lambda: 'mismatch')
    for line in pancan_fh: 
        line = line.strip().split('\t')
        pt = line[0]
        val = line[1]
        pancan_dict[pt] = val
    
    expr_fh = open(args.expr)
    expr_head = expr_fh.next().strip().split('\t')
    
    headero = [
               'transcript_id',
#                'gene_id',
               'pt_id',
               'can_val',
               'type_val',
               'pancan_val',
               'expr_val'
               ]
    print '\t'.join(headero)
    for line in expr_fh: 
        line = line.strip().split('\t')
        t_id = line[expr_head.index('transcript_id')]
        #g_id = line[expr_head.index('gene_id')]
        for x in xrange(LAST_PHENO_COL, len(line)):
            expr_val = line[x]
            pt = expr_head[x]
            peak_val = canc_dict[pt]
            type_val = type_dict[pt]
            pancan_val = pancan_dict[pt]
            lineo = t_id, pt, peak_val, type_val, pancan_val,expr_val
            print '\t'.join(map(str, lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
