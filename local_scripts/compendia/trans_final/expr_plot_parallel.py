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


def worker(args):
    (expr_mat, out_dir, i) = args
    
    if i%500==0:
        logging.debug('On worker %d' %i)
    
    r_args = [
              'Rscript',
              '/mctp/projects/mitranscriptome/naming/all_type_expr_fpkm.R',
              expr_mat,
              out_dir
              ]
    
    subprocess.call(r_args)
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("mat_dir")
    parser.add_argument("out_dir")
    parser.add_argument("-p", dest = 'proc',
                    default = 4,
                    help = 'number of processors to use')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    expr_mats = glob.glob(os.path.join(args.mat_dir, '*'))
    
    
    pool = multiprocessing.Pool(int(args.proc))
    tasks = []
    i=0
    for expr_mat in expr_mats:
        i+=1
        r_args = (expr_mat, args.out_dir, i)
        tasks.append(r_args)
    pool.map(worker, tasks)
    pool.close()
    pool.join()
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
