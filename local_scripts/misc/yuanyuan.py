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


del_meta_col = 1
expr_meta_col = 13
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("del_file")
    parser.add_argument("expr_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    del_file = open(args.del_file)
    del_header = del_file.next().strip().split('\t')
    del_meta_header = del_header[:del_meta_col]
    del_vals_header = del_header[del_meta_col:]
    dlc_del_dict = collections.defaultdict(lambda: 'NA')
    erg_del_dict = collections.defaultdict(lambda: 'NA')
    for line in del_file: 
        line = line.strip().split('\t')
        line_meta = line[:del_meta_col]
        line_vals = line[del_meta_col:]
        gene_name = line[del_header.index("Gene Symbol")]
        if gene_name == 'DLC1':
            for x in xrange(len(del_vals_header)): 
                pt = del_vals_header[x]
                del_val = line_vals[x]
                dlc_del_dict[pt] = del_val
        if gene_name == 'ETV1':
            for x in xrange(len(del_vals_header)): 
                pt = del_vals_header[x]
                del_val = line_vals[x]
                erg_del_dict[pt] = del_val
    
    expr_file = open(args.expr_file)
    expr_header = expr_file.next().strip().split('\t')
    expr_meta_header = expr_header[:expr_meta_col]
    expr_vals_header = expr_header[expr_meta_col:]
    dlc_expr_dict = collections.defaultdict(lambda: 'NA')
    erg_expr_dict = collections.defaultdict(lambda: 'NA')
    for line in expr_file: 
        line = line.strip().split('\t')
        line_meta = line[:expr_meta_col]
        line_vals = line[expr_meta_col:]
        gene_name = line[expr_header.index("gene_name")]
        if gene_name == 'DLC1':
            for x in xrange(len(expr_vals_header)): 
                pt = expr_vals_header[x]
                expr_val = line_vals[x]
                dlc_expr_dict[pt] = expr_val
        if gene_name == 'ETV1':
            for x in xrange(len(expr_vals_header)): 
                pt = expr_vals_header[x]
                expr_val = line_vals[x]
                erg_expr_dict[pt] = expr_val
    
    headero = ['pt_id', 'dlc_del', 'dlc_expr', 'etv1_del', 'etv1_expr']
    print '\t'.join(headero)
    for pt in dlc_del_dict.iterkeys():
        dlc_del_val = dlc_del_dict[pt]
        dlc_expr_val = dlc_expr_dict[pt]
        erg_del_val = erg_del_dict[pt]
        erg_expr_val = erg_expr_dict[pt]
        if (dlc_del_val == 'NA') or (dlc_expr_val == 'NA') or (erg_del_val == 'NA') or (erg_expr_val == 'NA'):
            continue
        lineo = [pt, dlc_del_val, dlc_expr_val, erg_del_val, erg_expr_val]
        print '\t'.join(map(str, lineo))
        
    return 0

if __name__ == '__main__': 
    sys.exit(main())
