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
mut_meta_col = 1
expr_meta_col = 12

'''
will take a gene expression matrix, and a matrix containing copy number data 
and will output a list of ptID, expression of gene, copy number of gene for each patient

will also take a mutation file and report status of mutation for TP53

copy number file acquired from the UCSC cancer browser, "GISTIC" file. 
'''
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("del_file")
    parser.add_argument("expr_file")
    parser.add_argument("--mut", dest='mut')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    del_file = open(args.del_file)
    del_header = del_file.next().strip().split('\t')
    del_meta_header = del_header[:del_meta_col]
    del_vals_header = del_header[del_meta_col:]
    cdc20_del_dict = collections.defaultdict(lambda: 'NA')
    for line in del_file: 
        line = line.strip().split('\t')
        line_meta = line[:del_meta_col]
        line_vals = line[del_meta_col:]
        gene_name = line[del_header.index("Gene Symbol")]
        if gene_name == 'CDC20':
            for x in xrange(len(del_vals_header)):
                pt = del_vals_header[x]
                del_val = line_vals[x]
                cdc20_del_dict[pt] = del_val
    
    expr_file = open(args.expr_file)
    expr_header = expr_file.next().strip().split('\t')
    expr_meta_header = expr_header[:expr_meta_col]
    expr_vals_header = expr_header[expr_meta_col:]
    cdc20_expr_dict = collections.defaultdict(lambda: 'NA')
    for line in expr_file: 
        line = line.strip().split('\t')
        line_meta = line[:expr_meta_col]
        line_vals = line[expr_meta_col:]
        gene_name = line[expr_header.index("gene_name")]
        if gene_name == 'CDC20':
            for x in xrange(len(expr_vals_header)): 
                pt = expr_vals_header[x]
                expr_val = line_vals[x]
                cdc20_expr_dict[pt] = expr_val
    
    mut_file = open(args.mut)
    mut_header = mut_file.next().strip().split('\t')
    mut_meta_header = mut_header[:mut_meta_col]
    mut_vals_header = mut_header[mut_meta_col:]
    cdc20_mut_dict = collections.defaultdict(lambda: 'NA')
    for line in mut_file: 
        line = line.strip().split('\t')
        line_meta = line[:mut_meta_col]
        line_vals = line[mut_meta_col:]
        gene_name = line[mut_header.index("sample")]
        if gene_name == 'TP53':
            for x in xrange(len(mut_vals_header)): 
                pt = mut_vals_header[x]
                mut_val = line_vals[x]
                cdc20_mut_dict[pt] = mut_val
    
    
    headero = ['pt_id', 'cdc20_del', 'cdc20_expr', 'tp53_status']
    print '\t'.join(headero)
    for pt in cdc20_del_dict.iterkeys():
        cdc20_del_val = cdc20_del_dict[pt]
        cdc20_expr_val = cdc20_expr_dict[pt]
        cdc20_mut_val = cdc20_mut_dict[pt]
        if (cdc20_expr_val == 'NA') or (cdc20_mut_val == 'NA'):
            continue
        lineo = [pt, cdc20_del_val, cdc20_expr_val, cdc20_mut_val]
        print '\t'.join(map(str, lineo))
        
    return 0

if __name__ == '__main__': 
    sys.exit(main())
