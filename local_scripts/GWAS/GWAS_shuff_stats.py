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
from scipy import stats

'''
takes a directory with OR data for multiple flanks 
and provides a table of summary stats
'''


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    header = [
              'analysis',
              'mean_OR_gwas',
              'std_OR_gwas',
              'mean_OR_rand',
              'std_OR_rand',
              't_stat',
              'p_value'
              ]
    print '\t'.join(header)
    files = glob.glob(args.result_dir+'/*.txt')[:]
    for file in files:
        logging.debug(file)
        base = os.path.basename(file)[:-4]
        ORs = np.array([])
        ORs_rand = np.array([])
        fh = open(file)
        header = fh.next().strip().split('\t')
        for line in fh: 
            line = line.strip().split('\t')
            OR = float(line[header.index("OR_gwas")])
            OR_rand = float(line[header.index("OR_rand")])
            ORs = np.append(ORs, OR)
            ORs_rand = np.append(ORs_rand, OR_rand)
        mean = ORs.mean()
        mean_rand = ORs_rand.mean()
        std = np.std(ORs)
        std_rand = np.std(ORs_rand)
        t, pval = stats.ttest_rel(ORs, ORs_rand)
        lineo = [
                 base,
                 mean,
                 std,
                 mean_rand,
                 std_rand,
                 t,
                 pval
                 ]
        print '\t'.join(map(str, lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
