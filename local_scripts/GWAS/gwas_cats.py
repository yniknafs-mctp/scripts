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

'''
Bladder
Blood
Brain
Breast
Colorectal
Esophagus
Gallbladder
HN
Kidney
Liver
Lung
NA
Other
Ovarian
Pancreas
Prostate
Skin
Thyroid
Tissue
Uterus
'''


key_pairs = [
             ('Bladder','bladder'),
             ('Breast','breast'),
             ('Blood','aml'),
             ('Blood','cml'),
             ('Blood','mpn'),
             ('Brain','gbm'),
             ('Brain','lgg'),
             ('Colorectal','colorectal'),
             ('Skin','melanoma'),
             ('Uterus','uterine'),
             ('Kidney','kich'),
             ('Kidney','kirp'),
             ('Kidney','kirc'),
             ('HN','head_neck'),
             ('Liver','liver'),
             ('Lung','luad'),
             ('Lung','lusc'),
             ('Ovarian','ovarian'),
             ('Pancreas','pancreas'),
             ('Prostate','prostate'),
             ('Thyroid','thyroid')
             ]


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
    cscat_header.append('gwas_match')
    print '\t'.join(cscat_header)
    for line in cscat_fh: 
        can_type = []
        line = line.strip().split('\t')
        t_id = line[cscat_header.index('transcript_id')]
        can_type.append(line[cscat_header.index('func_type')])
        can_type.append(line[cscat_header.index('ssea_type.ctnt.up')])
        can_type.append(line[cscat_header.index('ssea_type.ctnt.dn')])
        can_type.append(line[cscat_header.index('ssea_type.cn.up')])
        can_type.append(line[cscat_header.index('ssea_type.cn.dn')])
        
        gwas_info = gwas_dict[t_id]
        gwas_cancer = gwas_type_dict[t_id]
        gwas_match=0
        for (a, b) in key_pairs:
            if gwas_cancer == a and b in can_type:
                gwas_match=1
                    

        for item in gwas_info:
            line.append(item)
        line.append(gwas_match)
        print '\t'.join(map(str, line))
            
        
        
        
        
        
    return 0

if __name__ == '__main__': 
    sys.exit(main())
