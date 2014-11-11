'''
Created on Jul 23, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess
import collections
import datetime

'''
Takes a file listing all of the tissue types and the locations
of the SMX files to generate tissue expression stats for each tissue type
'''

TRANS_META_SCRIPT = '/mctp/users/yniknafs/scripts/workspace_laptop/misc_scripts/local_scripts/rna_seq/trans_meta_expr_one_type.py'
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_dir")
    parser.add_argument("smx_files")
    parser.add_argument("smx_dir")
    parser.add_argument("metadata")
    parser.add_argument("out_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    TMP_DIR = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
    
    os.mkdir(TMP_DIR)
    
    tissue_dict = collections.defaultdict(lambda: [])
    
    for line in open(args.smx_files):
        tissue, smx = line.strip().split('\t')
        tissue_dict[tissue].append(smx)
    
    
    i=0
    for tissue, smxs in tissue_dict.iteritems():
        logging.debug('Calculating expression stats for %s' % tissue)
        i+=1
        LIB_TMP = os.path.join(TMP_DIR, 'lib_' + tissue + '.txt')
        with open(LIB_TMP, 'w') as f:
            print >>f, 'library_id'
            for smx in smxs: 
                smx_file = open(os.path.join(args.smx_dir, smx))
                smx_file.next()
                smx_file.next()
                for line in smx_file: 
                    line = line.strip().split('\t')
                    if line[1] == '1': print >>f, line[0]
        
        if i == 1: OUT_TMP1 = args.metadata 
        else: OUT_TMP1 = os.path.join(TMP_DIR, str(i-1) + '_tmp.txt')
        OUT_TMP2 = os.path.join(TMP_DIR, str(i) + '_tmp.txt')
        
        trans_meta_args = [
                           'python',
                           TRANS_META_SCRIPT,
                           args.matrix_dir,
                           tissue, 
                           LIB_TMP,
                           OUT_TMP1
                           ]
        with open(OUT_TMP2, 'w') as f:
            subprocess.call(trans_meta_args, stdout=f)
        
    os.rename(OUT_TMP2, os.path.abspath(args.out_file))
    rm_args = ['rm','-rf', TMP_DIR]
    subprocess.call(rm_args)
    logging.debug('Done')
        
        
        
    return 0 
        
        

if __name__ == '__main__': 
    sys.exit(main())
