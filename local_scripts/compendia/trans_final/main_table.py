'''
Created on Jun 2, 2014
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
    
    headero = [
               'tissue',
               'tot',
               'clat',
               'cons',
               'uce',
               'tucp',
               'gwas'
               ]
    
    file_h = open(args.meta)
    file_header = file_h.next().strip().split('\t')
    tissue_dict = collections.defaultdict(lambda: [])
    cons_dict = collections.defaultdict(lambda: [])
    uce_dict = collections.defaultdict(lambda: [])
    clat_dict = collections.defaultdict(lambda: [])
    tucp_dict = collections.defaultdict(lambda: [])
    gwas_dict = collections.defaultdict(lambda: [])
    print '\t'.join(headero)
    for line in file_h:
        line = line.strip().split('\t')
        gene_id = line[file_header.index('gene_id')]
        tissue = line[file_header.index('func_type')]
        cons = line[file_header.index('cons')]
        uce = line[file_header.index('uce')]
        clat = line[file_header.index('func_cat')]
        tcat = line[file_header.index('tcat')]
        gwas = line[file_header.index('gwas_snp_dist')]
        tissue_dict[tissue].append(gene_id)
        if cons == 'TRUE':
            cons_dict[tissue].append(gene_id)
        if uce == 'TRUE': 
            uce_dict[tissue].append(gene_id)
        if clat == 'Cancer and Lineage Association':
            clat_dict[tissue].append(gene_id)
        if tcat == 'tucp':
            tucp_dict[tissue].append(gene_id)
        if gwas == '0':
            gwas_dict[tissue].append(gene_id)
    for tissue in tissue_dict.iterkeys():
        genes = len(set(tissue_dict[tissue]))
        cons = len(set(cons_dict[tissue]))
        uce = len(set(uce_dict[tissue]))
        clat = len(set(clat_dict[tissue]))
        tucp = len(set(tucp_dict[tissue]))
        gwas = len(set(gwas_dict[tissue]))
        lineo = [tissue, genes, clat, cons, uce, tucp, gwas]
        print '\t'.join(map(str,lineo))
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
