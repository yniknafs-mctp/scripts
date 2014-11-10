'''
Created on Jan 18, 2014
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
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    cnv_script = '/home/yniknafs/workspace/local_scripts/CNV/gistic_proc/gistic_main_set_maker.py' 
    cnv_file = '/mctp/users/yniknafs/projects/CNV/ca_types.tsv'
   
    out_path = 'CNV_processed'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    
    for line in open(cnv_file):
        line = line.strip().split('\t')
        lesion_file = line[0]
        matrix_file = line[1]
        logging.debug('lesion file: %s' % lesion_file)
        logging.debug('matrix file: %s' % matrix_file)
        ca_type = line[2]
        logging.debug(ca_type)
        result_dir = os.path.join(out_path, ca_type)
        if not os.path.exists(result_dir): 
            os.mkdir(result_dir)
        
        logging.info('Running CNV processing for  %s' % ca_type)        
        
        args_cnv = ['python', cnv_script,
                      lesion_file,
                      matrix_file,
                      ca_type,
                      result_dir]
        
        for arg in args_cnv:
            if arg == 'python ': 
                continue
            if not os.path.exists(arg):
                logging.debug('%s does not exist' % arg)
        
        
        
        subprocess.call(args_cnv)
        

        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
