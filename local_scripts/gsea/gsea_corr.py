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
GSEA_PATH = '/mctp/users/yniknafs/sw/gsea2-2.1.0.jar'
RAM = '2000m'
SET_MAX = '1000'
SET_MIN = '10'
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


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("tid")  #transcript that is being correlated
    parser.add_argument("trans_meta")
    parser.add_argument("lib_meta")
    parser.add_argument("tissue_type")
    parser.add_argument("matrix_dir")
    parser.add_argument("--gmt", dest = "gmt")
    parser.add_argument("-o", dest = "out", default = None)
    parser.add_argument("--expr_mat_dir", dest = "expr_mat_dir", default = None)
    parser.add_argument("--perms", dest = "perms")
    parser.add_argument("--job_id", dest='jobid', default = '1')
    
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
        if tid == args.tid:
            ref_gene_list.append(tid)
        cat = line[header.index('category')]
        ref_name = line[header.index('ref_gene_name')].upper()
        if ((cat == 'same_strand') or (cat == 'read_through')):
            ref_name_dict[tid] = ref_name
            ref_gene_list.append(tid)
    
    TMP_DIR = os.path.join(args.out, 'GSEA_CORR_TMP' + args.jobid)
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
            if t_id ==args.tid:
                print >>f, '\t'.join(line)
    
    tissue_type = tissue_translate(args.tissue_type)
    if tissue_type == 'NA':
        logging.error('Not a valid tissue type: %s' % args.tissue_type)
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
        
    logging.debug("Acquiring expression data")
    TMP_EXPR_MAT = os.path.join(TMP_DIR, 'expr.txt')
    with open(TMP_EXPR_MAT, 'w') as f: 
        subprocess.call(bcm_to_tsv_args, stdout=f)
    
    expr_mat = os.path.join(args.expr_mat_dir, args.tissue_type + '.expr.txt')
    expr_mat_fh = open(expr_mat)
    expr_mat_header = expr_mat_fh.next()
    tmp_expr_mat_header = open(TMP_EXPR_MAT).next()
    if expr_mat_header != tmp_expr_mat_header:
        logging.error("Expression Matrices Do Not Match")
        return 1
    
    with open(TMP_EXPR_MAT, 'a') as f: 
        for line in expr_mat_fh:
            print >> f, line.strip()
    
    
    #read expression data to grab transcript of interest
    fileh = open(TMP_EXPR_MAT)
    fileh.next()
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0]
        if tid == args.tid:
            a = np.array(line[1:]).astype(np.float)
    #go back through file and determine correlations
    logging.info("Correlating...")
    fileh = open(TMP_EXPR_MAT)
    fileh.next()
    cors = []
    i = 0
    for line in fileh:
        i+=1
#         if i>50000:
#             continue
        if i%1000==0:
            logging.debug("finished %d" % i)
        line = line.strip().split('\t')
        tid = line[0]
        if tid == args.tid:
            continue
        ref_name = ref_name_dict[tid]
        b = np.array(line[1:]).astype(np.float)
        r,p = spearmanr(a, b)
        if np.isnan(r): 
            continue 
        for gene in ref_name.split(','):
            cors.append([gene, r])
         
    #go back through file and determine correlations
    logging.info("Sorting...")
    cor = sorted(cors, key=itemgetter(1), reverse=True)
    cor_r = sorted(cors, key=itemgetter(1))
    max_list = set()
    min_list = set()
    corr_final_pos = []
    #pick best gene isoform
    logging.info("Picking best isoform...")
    for gene, r in cor:
        if r<0:
            continue
        if gene in max_list:
            continue
        max_list.add(gene)
        corr_final_pos.append([gene, r])
    corr_final_neg = []
    for gene, r in cor_r:
        if r>0:
            continue
        if gene in min_list:
            continue
        min_list.add(gene)
        corr_final_neg.append([gene, r])    
     
    corr_final = corr_final_pos + corr_final_neg[::-1]
     
    RNK_LIST = os.path.join(TMP_DIR, 'cors.rnk')
    with open(RNK_LIST, 'w') as f: 
        for pair in corr_final:
            print>>f, '\t'.join(map(str, pair))
     
    logging.debug("Running GSEA")
    gsea_params = [
                   'java', 
                   '-cp', 
                   GSEA_PATH,
                   '-Xmx' + RAM, 
                   'xtools.gsea.GseaPreranked',
                   '-gmx',
                   args.gmt,
                   '-collapse',
                   'false',
                   '-mode',
                   'Max_probe',
                   '-norm',
                   'meandiv',
                   '-nperm',
                   args.perms,
                   '-rnk',
                   RNK_LIST,
                   '-scoring_scheme', 
                   'weighted',
                   '-rpt_label',
                   os.path.basename(TMP_DIR),
                   '-include_only_symbols',
                   'true',
                   '-make_sets',
                   'true',
                   '-plot_top_x',
                   '0',
                   '-rnd_seed',
                   'timestamp',
                   '-set_max',
                   SET_MAX,
                   '-set_min',
                   SET_MIN,
                   '-zip_report',
                   'false',
                   '-out',
                   TMP_DIR,
                   '-gui',
                   'false'                   
                   ]
    subprocess.call(gsea_params)
     
    GSEA_RES_DIR = glob.glob(os.path.join(TMP_DIR, '*.GseaPreranked*'))[0]
    gsea_results = glob.glob(os.path.join(GSEA_RES_DIR, '*'))
    
    if args.out:
        OUT_DIR = args.out
    else: 
        logging.error('No output directory specified')
        return 1 
    
    out_pref = args.tid + '_' + args.tissue_type 
    
    for file in gsea_results: 
        base = os.path.basename(file)
        if base.startswith('gsea_report_for_na_pos') and base.endswith('xls'):
            os.rename(file, os.path.join(OUT_DIR, out_pref + '_pos.txt'))
        if base.startswith('gsea_report_for_na_neg') and base.endswith('xls'):
            os.rename(file, os.path.join(OUT_DIR, out_pref +'_neg.txt'))
#         elif file.endswith('edb'): 
#             shutil.rmtree(file)
#         else: 
#             os.remove(file)
    shutil.rmtree(GSEA_RES_DIR)
    for fileh in glob.glob(TMP_DIR+"/*"):
        os.remove(fileh)
    subprocess.call('rm -rf %s' % TMP_DIR, shell = True)
    JOB_RUNNING = os.path.join(OUT_DIR, out_pref + '.running')
    JOB_DONE = os.path.join(OUT_DIR, out_pref + '.done')
    if os.path.exists(JOB_RUNNING):
        os.rename(JOB_RUNNING, JOB_DONE)    
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

