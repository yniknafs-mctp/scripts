'''
Created on Apr 16, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
from random import shuffle


# EXP_LEVELS = [1,100]
# FOLD_LEVELS = [1,50]
# OUTLIER_LEVELS = [1,.5,.05]



EXP_LEVELS = [1,10,25,100]
FOLD_LEVELS = [1,2,5,7,10,25,50]
OUTLIER_LEVELS = [1,.5,.35 ,.30 ,.25, .20, .15, .1, .05]

DISP_LIMITS = [50,200]
DIRS = ['up','dn']
FOLD_CHANGE_FRAC = .5
REPS_FC = 500
TOT_FC = len(EXP_LEVELS)*(len(FOLD_LEVELS)-1)*len(OUTLIER_LEVELS)*REPS_FC
TOT_CTRL = TOT_FC/FOLD_CHANGE_FRAC - TOT_FC
REPS_CTRL = int(TOT_CTRL/(len(EXP_LEVELS)*len(OUTLIER_LEVELS)))
print TOT_FC
print TOT_CTRL



def count_transform(count, method, size):
    if method == 'uniform': 
        return count
    if method == 'poisson':
        return np.random.poisson(count,1)[0]
    if method == 'nb': 
        prob = float(size)/(size+count)
        return np.random.negative_binomial(size, prob, 1)[0]
    else: 
        print 'not compatible method'
        sys.exit()
        
def nb(count, size):
    prob = size/(size+count)
    return np.random.negative_binomial(1, prob, size)[0]
        
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_meta")
    parser.add_argument("expr")
    parser.add_argument("method")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    TRANSFORM = args.method
    
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    full_mat = os.path.join(output_dir, 'full_mat.txt')
    sim_mat = os.path.join(output_dir, 'sim_mat.txt')
    sim_meta = os.path.join(output_dir, 'sim_meta.txt')
    
    with open(full_mat, 'w') as f_full: 
        expr_fh = open(args.expr)
        expr_header = expr_fh.next().strip().split('\t')
        print >>f_full, '\t'.join(expr_header)
#         logging.debug('Scanning origninal file for copy')
#         for line in expr_fh: 
#             line = line.strip()
#             print >>f_full, line
            
        sample_meta_fh = open(args.sample_meta)
        sample_meta_header = sample_meta_fh.next().strip().split('\t')
        norms = []
        cans = []
        for line in sample_meta_fh: 
            line = line.strip().split('\t')
            lib_id = line[sample_meta_header.index('library_id')]
            prog = line[sample_meta_header.index('cancer_progression')]
            if prog == 'normal': 
                norms.append(lib_id)
            else: 
                cans.append(lib_id)
        
        
        meta_header = [
                       'transcript_id',
                       'norm_expression',
                       'can_expression',
                       'fold_change',
                       'outlier_frac',
                       'direction',
                       'size'
                       ]

        i = 0
        logging.debug('Simulating data')
        with open(sim_mat, 'w') as f_sim:
            print >>f_sim, '\t'.join(expr_header)
            with open(sim_meta, 'w') as f_meta:
                print >>f_meta, '\t'.join(meta_header)
                for exp in EXP_LEVELS: 
                    shuffle(cans)
                    for fold in FOLD_LEVELS:
                        for frac in OUTLIER_LEVELS: 
                            for dir in DIRS:
                                if fold==1:
                                    REPS = int(REPS_CTRL)
                                else: 
                                    REPS = int(REPS_FC)
                                for x in xrange(REPS):
                                    i+=1
                                    if i%500==0:
                                        logging.debug('Finished %d simulations' %(i))
                                    num_can = int(frac*len(cans))
                                    exp_dict = {}
                                    shuffle(cans)
                                    seed = exp*fold
                                    name = 'TS' + str(i)
                                    size = np.random.uniform(DISP_LIMITS[0], DISP_LIMITS[1],1)[0]
                                    if dir == 'up': 
                                        can_exp = seed
                                        norm_exp = exp
                                    if dir == 'dn':
                                        can_exp = exp
                                        norm_exp = seed
                                    for pt in cans[:num_can]: 
                                        exp_dict[pt] = count_transform(can_exp, TRANSFORM, size) 
                                    for pt in cans[num_can:]: 
                                        exp_dict[pt] = count_transform(norm_exp, TRANSFORM, size)
                                    for pt in norms: 
                                        exp_dict[pt] = count_transform(norm_exp, TRANSFORM, size)
                                    lineo = [name]
                                    for pt in expr_header[1:]:
                                        expr = exp_dict[pt]
                                        lineo.append(expr)
                                    lineo_meta = [name, norm_exp, can_exp, fold, frac, dir, size]
                                    print >>f_sim, '\t'.join(map(str,lineo))
                                    print >>f_full, '\t'.join(map(str,lineo))
                                    print >>f_meta, '\t'.join(map(str, lineo_meta))
                
            
    return 0

if __name__ == '__main__': 
    sys.exit(main())
