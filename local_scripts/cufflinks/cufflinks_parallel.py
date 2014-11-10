'''
Created on May 10, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess
import multiprocessing
import glob

COMPUTE_GENE_EXPR_SCRIPT = '/mctp/users/yniknafs/scripts/workspace/assemblyline/assemblyline/pipeline/compute_gene_expression.py'

def worker(args):
    (sample,tot, i) = args 
    logging.info('Running cufflinks for %s; worker %d/%d' % (os.path.basename(sample), i, tot))
    sh_args = [
               'sh',
               os.path.join(sample, 'cufflinks.sh')
              ]
    with open(os.path.join(sample, 'log.txt'), 'w') as f:
        subprocess.call(sh_args, stderr=f)
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml")
    parser.add_argument("library_table")
    parser.add_argument("gtf_file")
    parser.add_argument("out_dir")
    parser.add_argument("-p", dest = 'proc',
                    default = 4,
                    help = 'number of processors to use')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #run the script to make the PBS scripts
    script_gen_args = [
                       'python',
                       COMPUTE_GENE_EXPR_SCRIPT,
                       '-o',
                       args.out_dir, 
                       args.config_xml,
                       args.library_table,
                       args.gtf_file,
                       '--mode',
                       'cufflinks' 
                       ]
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    
    
    subprocess.call(script_gen_args)
    
    sample_dirs = glob.glob(os.path.join(args.out_dir, '*'))
    
    #multiprocess run cufflinks
    pool = multiprocessing.Pool(int(args.proc))
    tasks = []
    i=0
    for sample in sample_dirs:
        i+=1
        cuff_args = (sample, len(sample_dirs), i)
        tasks.append(cuff_args)
    pool.map(worker, tasks)
    pool.close()
    pool.join()
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
