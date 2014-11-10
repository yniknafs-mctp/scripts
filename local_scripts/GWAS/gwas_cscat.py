'''
Created on Apr 15, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np


'''
merge cscat and gwas info
return cscats that overlap snps of the correct type
'''


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("cscat_meta")
    parser.add_argument("gwas_meta")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    gwas_fh = open(args.gwas_meta)
    gwas_header = gwas_fh.next().strip().split('\t')
    
    gwas_dict = {}
    gwas_type_dict = {}
    logging.debug('Obtaining gwas data')
    for line in gwas_fh:
        gwas_info = [] 
        line = line.strip().split('\t')
        t_id = line[gwas_header.index('transcript_id')]
        gwas_snp = line[gwas_header.index('gwas_snp_id')]
        gwas_dist = line[gwas_header.index('gwas_snp_dist')]
        gwas_exon = line[gwas_header.index('gwas_exonic_overlap')]
        gwas_cancer = line[gwas_header.index('gwas_cancer')]
        gwas_info.append(gwas_snp)
        gwas_info.append(gwas_dist)
        gwas_info.append(gwas_exon)
        gwas_info.append(gwas_cancer)
        gwas_dict[t_id] = gwas_info
        gwas_type_dict[t_id] = gwas_cancer
    
    logging.debug('Searching for compendia trnascripts that overlap')
    cscat_fh = open(args.cscat_meta)
    cscat_header = cscat_fh.next().strip().split('\t')
    cscat_header.append('gwas_snp_id')
    cscat_header.append('gwas_snp_dist')
    cscat_header.append('gwas_exonic_overlap')
    cscat_header.append('gwas_cancer')
    cscat_header.append('gwas_cscat')
    print '\t'.join(cscat_header)
    for line in cscat_fh: 
        line = line.strip().split('\t')
        t_id = line[cscat_header.index('transcript_id')]
        cscat_up = line[cscat_header.index('ssea_type_up')]
        type_up = line[cscat_header.index('cancer_type_up')]
        gwas_info = gwas_dict[t_id]
        gwas_cancer = gwas_type_dict[t_id]
        gwas_cscat=0
        if cscat_up=='cscat':
            if gwas_cancer == 'Bladder' and type_up == 'bladder':
                gwas_cscat=1    
            if gwas_cancer == 'Breast' and type_up == 'breast':
                gwas_cscat=1
            if gwas_cancer == 'Kidney' and type_up == 'kich':
                gwas_cscat=1
            if gwas_cancer == 'Kidney' and type_up == 'kirc':
                gwas_cscat=1
            if gwas_cancer == 'Kidney' and type_up == 'kirp':
                gwas_cscat=1
            if gwas_cancer == 'HN' and type_up == 'head_neck':
                gwas_cscat=1
            if gwas_cancer == 'Liver' and type_up == 'liver':
                gwas_cscat=1
            if gwas_cancer == 'Lung' and type_up == 'luad':
                gwas_cscat=1
            if gwas_cancer == 'Lung' and type_up == 'lusc':
                gwas_cscat=1
            if gwas_cancer == 'Ovarian' and type_up == 'ovarian':
                gwas_cscat=1
            if gwas_cancer == 'Pancreas' and type_up == 'pancreas':
                gwas_cscat=1
            if gwas_cancer == 'Prostate' and type_up == 'prostate':
                gwas_cscat=1
            if gwas_cancer == 'Thyroid' and type_up == 'thyroid':
                gwas_cscat=1
        for item in gwas_info:
            line.append(item)
        line.append(gwas_cscat)
        print '\t'.join(map(str, line))
            
        
        
        
        
        
    return 0

if __name__ == '__main__': 
    sys.exit(main())
