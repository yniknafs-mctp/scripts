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
    
    shuff_script = '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_intergenic_shuffle_parallel.py' 
    gwas_file = '/mctp/projects/mitranscriptome/gwas/bed_files/gwas.sorted.bed'
#     test_script = '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/test.py'
    
    
    '''
    flanks = [(100, '0.1'),
              (300, '0.3'),
              (500, '0.5'),
              (1000, '1'),
              (2000, '2'),
              (3000, '3'),
              (4000, '4'), 
              (5000, '5'),
              (7500, '7.5'),
              (10000, '10'),
              (20000, '20'),
              (30000, '30'),
              (50000, '50'),
              (75000, '75'),
              (100000, '100'),
              (125000, '125'),
              (150000, '150'),
              (200000, '200'),
              (500000, '500'),
              (750000, '750'),
              (1000000, '1000')]
    '''
    flanks = [(0, '0')]
    for flank, title in flanks:
        logging.info('Running flank size %d' % flank)
        flank_file =  title + 'kb_intergenic_shuffle.txt'
        
        args_shuff = ['python', shuff_script, 
                        '--shuffs', str(250),
                        '-p', str(60), 
                        '--gwas', gwas_file, 
                        '--flank', str(flank)]
        
        
        with open(flank_file, 'w') as fileh:
            subprocess.call(args_shuff, stdout=fileh)
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
