'''
Created on Jan 20, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import glob

'''
takes a directory with OR data for multiple flanks 
and provides a table of summary stats
'''


ORCOL = 3


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("flanks_dir")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    header = [
              'flank_size',
              'mean_OR',
              'std',
              'median_OR',
              'max_OR',
              'min_OR',
              '5th_tile',
              '95th_tile'
              ]
    print '\t'.join(header)
    for file in glob.glob(args.flanks_dir+'*.txt'):
        base = os.path.basename(file)
        flank_size = int(float(base.split('kb')[0])*1000)
        ORs = np.array([])
        for line in open(file): 
            line = line.strip().split('\t')
            OR = float(line[ORCOL])
            ORs = np.append(ORs, OR)
        mean = ORs.mean()
        median = np.median(ORs)
        std = np.std(ORs)
        max = ORs.max()
        min = ORs.min()
        five_tile = np.percentile(ORs, 5)
        ninetyfive_tile = np.percentile(ORs, 95)
        lineo = [
                 flank_size,
                 mean,
                 std,
                 median,
                 max,
                 min,
                 five_tile,
                 ninetyfive_tile
                 ]
        print '\t'.join(map(str, lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
