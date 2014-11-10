'''
Created on Jan 17, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np

'''
Takes snps bed file (merged affy and illumina snp arrays) 
and chooses a number of random snps, totalling the number 
of GWAS snps
'''  




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_bed")
    parser.add_argument("snps_bed")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Counting number of GWAS SNPS')
    num_gwas = len(open(args.gwas_bed).readlines())
    logging.info('Counting number of available SNPS')
    avail_snps = open(args.snps_bed).readlines()
    num_snps = len(avail_snps)
    logging.info('Shuffling snps')
    nums = np.arange(num_snps)
    np.random.shuffle(nums)
    rand_snp_index = nums[1:num_gwas]
    for index in rand_snp_index: 
        print avail_snps[index].strip()
    
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
