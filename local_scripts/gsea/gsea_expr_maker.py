'''
Created on Aug 26, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
from scipy import stats
from scipy.stats import spearmanr
import collections
import subprocess
from operator import itemgetter
import glob
import shutil

BCM_TO_TSV_PATH = '/mctp/users/yniknafs/scripts/workspace/local_scripts/rna_seq/bcm_to_tsv.py'

def tissue_translate(tissue):
    if tissue == 'breast':
        return ['BRCA_cancer', 'BRCA_met']
    if tissue == 'cervical':
        return ['CESC_cancer']
    if tissue == 'colorectal':
        return ['COADREAD_cancer', 'COADREAD_met']
    if tissue == 'gbm':
        return ['GMB_cancer']
    if tissue == 'head_neck':
        return ['HNSC_cancer']
    if tissue == 'liver':
        return ['LIHC_cancer']
    if tissue == 'lgg':
        return ['LGG_cancer']
    if tissue == 'luad':
        return ['LUAD_cancer', 'LUAD_met']
    if tissue == 'lusc':
        return ['LUSC_cancer']
    if tissue == 'pancreatic':
        return ['PAAD_cancer']
    if tissue == 'prostate':
        return ['PRAD_cancer', 'PRAD_met']
    if tissue == 'kich':
        return ['KICH_cancer']
    if tissue == 'kirc':
        return ['KIRC_cancer']
    if tissue == 'kirp':
        return ['KIRP_cancer']
    if tissue == 'stomach':
        return ['STAD_cancer']
    if tissue == 'uterine':
        return ['UCEC_cancer']
    if tissue == 'melanoma':
        return ['SKCM_cancer', 'SKCM_met']
    if tissue == 'thyroid':
        return ['THCA_cancer',]
    if tissue == 'ovarian':
        return ['OV_cancer',]
    else:
        return 'NA'

TISSUES = ['breast','cervical','colorectal',
           'gbm','head_neck','liver',
           'lgg','luad','lusc','pancreatic',
           'prostate','kich','kirc',
           'kirp','stomach','uterine',
           'thyroid', 'melanoma', 'ovarian']
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("trans_meta")
    parser.add_argument("lib_meta")
    parser.add_argument("matrix_dir")
    parser.add_argument("-o", dest = "out", default = None)
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    if not os.path.exists(args.out):
        os.mkdir(args.out)
    
    fileh = open(args.trans_meta)
    header = fileh.next().strip().split('\t')
    ref_name_dict = collections.defaultdict(lambda: 'NA')
    ref_gene_list = []
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[0]
        cat = line[header.index('category')]
        ref_name = line[header.index('ref_gene_name')].upper()
        if ((cat == 'same_strand') or (cat == 'read_through')):
            ref_name_dict[tid] = ref_name
            ref_gene_list.append(tid)
    
    TMP_DIR = os.path.join(args.out, 'GSEA_EXPR_TMP')
    if not os.path.exists(TMP_DIR):
        os.makedirs(TMP_DIR)
    TMP_TRANS_META = os.path.join(TMP_DIR, 'trans_meta.txt')
    TMP_LIB_META = os.path.join(TMP_DIR, 'lib_meta.txt')
    
    ref_gene_list = set(ref_gene_list)
    
    with open(TMP_TRANS_META, 'w') as f: 
        trans_meta_fh = open(args.trans_meta)
        trans_meta_header = trans_meta_fh.next().strip()
        print >>f, trans_meta_header
        for line in trans_meta_fh: 
            line = line.strip().split('\t')
            t_id = line[trans_meta_header.index('transcript_id')]
            if t_id in ref_gene_list:
                print >>f, '\t'.join(line)
    
    for tissue in TISSUES:             
        tissue_type = tissue_translate(tissue)
        if tissue_type == 'NA':
            logging.error('Not a valid tissue type: %s' % tissue)
            return 1
        with open(TMP_LIB_META, 'w') as f: 
            lib_meta_fh = open(args.lib_meta)
            lib_meta_header = lib_meta_fh.next().strip().split('\t')
            print >>f, '\t'.join(lib_meta_header)
            for line in lib_meta_fh: 
                line = line.replace('\n', '').split('\t')
                can_type = line[lib_meta_header.index('all_sample_plot_id')]
                if can_type in tissue_type:
                    print >>f, '\t'.join(line)
            
        bcm_to_tsv_args = [
                               'python',
                               BCM_TO_TSV_PATH, 
                               '--fpkm',
                               '--colmeta',
                               TMP_LIB_META,
                               '--rowmeta',
                               TMP_TRANS_META,
                               args.matrix_dir
                               ]
            
        logging.debug("Acquiring expression data for %s" % tissue)
        EXPR_MAT = os.path.join(args.out, '%s.expr.txt' % tissue)
        with open(EXPR_MAT, 'w') as f: 
            subprocess.call(bcm_to_tsv_args, stdout=f)
    
    for fileh in glob.glob(TMP_DIR+"/*"):
        os.remove(fileh)
    subprocess.call('rm -rf %s' % TMP_DIR, shell = True)
        
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

