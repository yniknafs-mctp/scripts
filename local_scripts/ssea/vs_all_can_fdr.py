'''
Created on Feb 12, 2014
@author yniknafs
'''

import os
import sys
import argparse
import logging
import numpy as np
import json
import ssea.lib.cfisher as fisher

'''
Takes file containing the path for all set results 
to be used in the set stats report 
'''

header = [
          'ss_name',
          'fdr_cutoff',
          'sensitivity',
          'specificity',
          'ca_sig',
          'ca_not_sig',
          'all_sig',
          'all_not_sig',
          'true_pos',
          'false_pos',
          'false_neg',
          'true_neg',
          'odds_ratio',
          'fisher_p'
          ]

fdr_cutoffs = [1e-2, 1e-3, 1e-4, 1e-5]

def fisher_test(a, b, anull, bnull):
    
    # subset both sets
    null = anull.intersection(bnull)
    # find overlap
    tp = len(a.intersection(b))
    fp = len(a.difference(b))
    fn = len(b.difference(a))
    tn = len(null)
    # fisher exact test (one-sided hypothesis that LE is enricheD)
    fisher_p_value = fisher.pvalue(tp, fp, fn, tn).right_tail
    if (fp == 0) or (fn == 0):
        odds_ratio = float('inf')
    else:
        odds_ratio = (tp * tn) / float(fp * fn)
    return [tp, fp, fn, tn, odds_ratio, fisher_p_value]

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("sets_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
        
    print >> sys.stdout, '\t'.join(header)
    # read through SS file to obtain t_ids for both sets being compared
    tot = len(open(args.sets_file).readlines())
    sets_file = open(args.sets_file)
    j = 0
    for line in sets_file:
        j+=1
        line = line.strip().split('\t')
        ss_ca, ss_all = line
        
        ss_ca = ss_ca.strip()
        ss_all = ss_all.strip()
        #read sample_set json to obtain sample set name
        ca_sample_set = json.load(open(os.path.join(ss_ca, 'sample_set.json')))
        ca_ss_name = ca_sample_set['name']
        all_sample_set = json.load(open(os.path.join(ss_all, 'sample_set.json')))
        all_ss_name = all_sample_set['name']
        logging.info('Processing sample set %d/%d: %s and %s' % (j, tot, ca_ss_name, all_ss_name))
        #read through ca vs normal result json to obtain result stats
        ca_t_ids = []
        ca_fdrs = []
        ca_results = open(os.path.join(ss_ca, 'results.json'))
        ca_results_read = ca_results.readlines()
        logging.debug("Reading cancer vs normal file")
        i = 0
        for x in ca_results_read:
            i+=1
            if (i%25000) == 0: 
                logging.debug('Finished %d/%d transcripts' % (i, len(ca_results_read)))
            result = json.loads(x.strip())
            t_id = result['t_id']
            ca_t_ids.append(t_id)
            ca_fdrs.append(result['ss_fdr_q_value'])
        
        #read through ca vs all normal result json to obtain result stats
        all_t_ids = []
        all_fdrs = []
        all_results = open(os.path.join(ss_all, 'results.json'))
        all_results_read = all_results.readlines()
        logging.debug("Reading cancer vs all normal file")
        i = 0
        for x in all_results_read:
            i+=1
            if (i%25000) == 0: 
                logging.debug('Finished %d/%d transcripts' % (i, len(all_results_read)))
            result = json.loads(x.strip())            
            t_id = result['t_id']
            all_t_ids.append(t_id)
            all_fdrs.append(result['ss_fdr_q_value'])
            
        ca_fdrs = np.array(ca_fdrs)
        ca_t_ids = np.array(ca_t_ids)
        all_fdrs = np.array(all_fdrs)
        all_t_ids = np.array(all_t_ids)
        for cutoff in fdr_cutoffs: 
            ca_lt = ca_fdrs<cutoff
            ca_lt_ts = set(ca_t_ids[ca_lt])
            ca_gt = ca_fdrs>=cutoff
            ca_gt_ts = set(ca_t_ids[ca_gt])
            all_lt = all_fdrs<cutoff
            all_lt_ts = set(all_t_ids[all_lt])
            all_gt = all_fdrs>=cutoff
            all_gt_ts = set(all_t_ids[all_gt])
            fisher = fisher_test(all_lt_ts, ca_lt_ts, all_gt_ts, ca_gt_ts)
            tp = fisher[0]
            tn = fisher[3]
            sens = float(tp)/sum(ca_lt)
            spec = float(tn)/sum(ca_gt)
            lineo = [
                     ca_ss_name,
                     cutoff,
                     sens, 
                     spec,
                     sum(ca_lt), 
                     sum(ca_gt),
                     sum(all_lt),
                     sum(all_gt)
                     ]
            lineo = lineo + fisher
            logging.info('\t'.join(map(str, lineo)))
            print >> sys.stdout, '\t'.join(map(str, lineo))
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
