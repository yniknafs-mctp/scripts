'''
Created on Apr 6, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections
import glob

'''
concatenate many SMX files with an expression file 
'''


LAST_PHENO_COL = 1

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("smx_dir")
    parser.add_argument("expr")
    parser.add_argument("out_dir")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        
    
    big_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: 'mismatch'))
    files = glob.glob(os.path.join(args.smx_dir,'*.smx'))
    for smx in files:
        fh = open(smx)
        fh.next()
        fh.next()
        for line in fh:
            name = os.path.basename(smx) 
            line = line.strip().split('\t')
            pt = line[0]
            val = line[1]
            big_dict[name[5:-4]][pt] = val

    logging.debug(files)
    
    expr_fh = open(args.expr)
    expr_head = expr_fh.next().strip().split('\t')
    
    headero = [
               'transcript_id',
               'pt_id',
               'expr_val'
               ]
    for key in big_dict.iterkeys():
        headero.append(key)
    
    print '\t'.join(headero)
    for line in expr_fh: 
        line = line.strip().split('\t')
        t_id = line[expr_head.index('transcript_id')]
        fileo = os.path.join(args.out_dir, (t_id + '_mat.txt'))
        with open(fileo, 'w') as f:
            print >>f, '\t'.join(headero)
            for x in xrange(LAST_PHENO_COL, len(line)):
                expr_val = line[x]
                pt = expr_head[x]
                lineo = [t_id, pt, expr_val]
                for key in big_dict.iterkeys():
                    val = big_dict[key][pt]
                    lineo.append(val)
                print >>f, '\t'.join(map(str, lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
