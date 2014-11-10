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
    
    gwas_script = '/mctp/users/yniknafs/scripts/workspace/local_scripts/GWAS/GWAS_intergenic_real.py' 
    gwas_file = '/mctp/projects/mitranscriptome/gwas/bed_files/gwas.sorted.bed'
    
    
    flanks = [(0, '0'),
              (100, '0.1'),
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
    for flank, title in flanks:
        logging.info('Running flank size %d' % flank)
        
        args_gwas = ['python', gwas_script, 
                        '--gwas', gwas_file, 
                        '--flank', str(flank)]
        
        
        subprocess.call(args_gwas)
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
