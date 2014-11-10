'''
Created on Mar 31, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections
from operator import itemgetter

'''
takes an aggregated cancer_vs_normal file and for any transcript 
with a significant result in pancan_vs_norm will add a column that reports
a comma delimited list of all categories for which that transcript was significant 
'''

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("cvn_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    
    fileh = open(args.cvn_file)
    header = fileh.next().strip().split('\t')
    pancans = []
    ups = []
    dns = []
    scores_up_dict = collections.defaultdict(lambda: [['NA','NA','NA']])
    scores_dn_dict = collections.defaultdict(lambda: [['NA','NA','NA']])
    pancan_es_dict = {}
    pancan_fdr_dict = {}
    for line in fileh: 
        line = line.strip().split('\t')
        ca_type = line[header.index('cancer_type')]
        es = float(line[header.index('es')])
        frac = float(line[header.index('frac')])
        fdr = float(line[header.index('fdr')])
        t_id = line[header.index('transcript_id')]
        stats = [ca_type,frac, es]
        if es >0:
            scores_up_dict[t_id].append(stats)
        elif es <=0:
            scores_dn_dict[t_id].append(stats)
        if ca_type == 'pancancer':
            pancans.append(t_id)
            pancan_es_dict[t_id] = es
            pancan_fdr_dict[t_id] = fdr
            if es >0: 
                ups.append(t_id)
            elif es <=0: 
                dns.append(t_id)
        
        
    pancans = set(pancans)
    ups = set(ups)
    dns = set(dns)
    #print header 
    headero = [
               'transcript_id',
#                'sig_types_up',
#                'sig_types_dn',
               'pancan_direction',
               'pancan_es',
               'pancan_fdr'
               ]
    print '\t'.join(headero)
    for t_id in pancans:
        scores_up = scores_up_dict[t_id]
        scores_dn = scores_dn_dict[t_id]
        if len(scores_up) > 1: 
            scores_up.pop(0)
        if len(scores_dn) > 1: 
            scores_dn.pop(0)
        pancan_score = pancan_es_dict[t_id]
        pancan_fdr = pancan_fdr_dict[t_id]
        scores_up_sort = sorted(scores_up, key=itemgetter(1), reverse=True)
        scores_up_str = ','.join(':'.join(map(str, l)) for l in scores_up_sort)
        scores_dn_sort = sorted(scores_dn, key=itemgetter(1))
        scores_dn_str = ','.join(':'.join(map(str, l)) for l in scores_dn_sort)
        if t_id in ups:
            direction = 'up'
        elif t_id in dns:
            direction = 'dn'
#         lineo = [t_id, scores_up_str, scores_dn_str, direction,pancan_score]
        lineo = [t_id, direction,pancan_score, pancan_fdr]
        print '\t'.join(map(str, lineo))

    return 0

if __name__ == '__main__': 
    sys.exit(main())
