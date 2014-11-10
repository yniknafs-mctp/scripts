'''
Created on May 11, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("manuscript_results")
    parser.add_argument("num_procs")
    parser.add_argument("transcript_meta")
    parser.add_argument("out_dir")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    tot = len(open(args.manuscript_results).readlines())
    
    i = 0 
    for dir in open(args.manuscript_results):
        i+=1
        dir = dir.strip()
        basename = os.path.basename(dir)
        outname = os.path.join(args.out_dir, basename)
        logging.debug('Processing %d/%d: %s' % (i, tot , basename))
        report_args = [
                   'python',
                   '/mctp/users/yniknafs/scripts/git/ssea/ssea/ssea/utils/ssea_report',
                   dir,
                   '/mctp/projects/ssea/isoform_count_matrix_v7/',
                   '-o',
                   outname,
                   '-p',
                    args.num_procs,
                    '--png',
                    '--exprpng',
                    '--rowmeta',
                    args.transcript_meta
                   ]
        
        subprocess.call(report_args)
        
     
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
