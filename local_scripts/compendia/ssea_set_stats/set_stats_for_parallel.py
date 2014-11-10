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
import collections
import glob

'''
Takes file containing the path for all set results 
to be used in the set stats report 
'''


header = [
             'ss_name',
             'ss_cat',
             'ss_type',
             'frac_cutoff',
             'es_sign',
             'fdr_cutoff',
             'metric', 
             'threshold',
             'category',
             'count'
            ]

fdr_cutoffs = [1e-3, 1e-5, 0]
prec_cutoffs = [0.95, 0.75, 0.5, .25]
fo_cutoffs = [0.01, 0.05, 0.1, .25]
frac_cutoffs = [.99, .975, .95, .9]
categories = ['protein_coding', 'lncrna', 'pseudogene', 'mixed_readthrough', 'tucp']


def ss_categorize(ss_dirname):
    ss_name = os.path.basename(ss_dirname)
    split = ss_name.split('_')
    if len(split) == 2: 
        cat = 'clinical'
        type = split[0]
    elif ss_name.startswith('cancer_type'):
        cat = 'cancer_type'
        type = '_'.join(split[2:])
    elif ss_name.startswith('cancer_versus_all_normals'):
        cat = 'cancer_versus_all_normals'
        type = '_'.join(split[4:])
    elif ss_name.startswith('cancer_versus_normal'):
        cat = 'cancer_versus_normal'
        type = '_'.join(split[3:])
    elif ss_name.startswith('cell_type'):
        cat = 'normal_cell_type'
        type = '_'.join(split[2:])
    elif ss_name.startswith('normal_cell_type'):
        cat = 'normal_cell_type'
        type = '_'.join(split[3:])
    elif ss_name.startswith('normal_type'):
        cat = 'normal_type'
        type = '_'.join(split[2:])
    elif 'fusion' in split:
        cat = 'fusion'
        type = ss_name.split('_tophat_fusion_')[0]
    elif split[1] == 'amp':
        cat = 'amplification'
        type = split[0]    
    elif split[1] == 'del':
        cat = 'deletion'
        type = split[0]
    elif split[2] == 'mutation':
        cat = 'mutation'
        type = split[0]
    else: 
        cat = 'clinical'
        type = split[0]
    return cat, type
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_json")
    parser.add_argument("set_file")
    parser.add_argument("results_dir")
    parser.add_argument("out_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    #read transcript metadata file to create dictionary for 
    logging.info('reading metadata')
    meta_file = open(args.metadata_json)
    meta_dict = {}
    k=0
    for line in meta_file:
        k+=1
        if k%10000==0:
            logging.debug('Parsed %d metadata results' %k) 
        line_json = json.loads(line.strip())
        meta_dict[line_json['_id']] = line_json['tcat']
    # parse results dir for list of sets
    sets = []
    for line in open(args.set_file):
        line = line.strip()
        sets.append(os.path.join(args.results_dir, line))
    tot = len(sets)
    with open(args.out_file, 'w') as f:
        i = 0
        for ss in sets:
            i+=1
            #read sample_set json to obtain sample set name
            sample_set = json.load(open(os.path.join(ss, 'sample_set.json')))
            ss_name = sample_set['name']
            logging.info('Processing sample set %d/%d: %s' % (i, tot, ss_name))
            counts_dict = collections.defaultdict(lambda: 0)
            ss_cat, ss_type = ss_categorize(ss)
            #read through result json to obtain result stats
            results = open(os.path.join(ss, 'results.json'))
            j=0
            for x in results:
                j+=1
                if j%10000 == 0: 
                    logging.debug('Parsed %d results: %s' % (j, ss_name))
                result = json.loads(x.strip())
                t_id = result['t_id']
                t_cat = meta_dict[t_id]
                fp = result['core_misses']
                tp = result['core_hits']
                tn = result['null_misses']
                frac = abs(result['ss_frac'])
                es = result['es']
                fdr = result['ss_fdr_q_value']
                if fp==0 and tn == 0:                 
                    continue
                #loop through and count stats for all different categories
                fo = float(fp)/(fp+tn)
                prec = float(tp)/(tp+fp)
                for frac_co in frac_cutoffs: 
                    for fdr_co in fdr_cutoffs: 
                        for cat in categories: 
                            for prec_co in prec_cutoffs: 
                                key = (ss_name, frac_co, 'neg', fdr_co, 'precision', prec_co, cat)
                                counts_dict[key] +=0
                                if es<0 and frac>frac_co and fdr<=fdr_co and t_cat==cat and prec>prec_co:
                                    key = (ss_name, frac_co, 'neg', fdr_co, 'precision', prec_co, cat)
                                    counts_dict[key] +=1
                                key = (ss_name, frac_co, 'pos', fdr_co, 'precision', prec_co, cat)
                                counts_dict[key] +=0
                                if es>=0 and frac>frac_co and fdr<=fdr_co and t_cat==cat and prec>prec_co:
                                    key = (ss_name, frac_co, 'pos', fdr_co, 'precision', prec_co, cat)
                                    counts_dict[key] +=1
                            for fo_co in fo_cutoffs:
                                key = (ss_name, frac_co, 'neg', fdr_co, 'fall_out', fo_co, cat)
                                counts_dict[key] +=0
                                if es<0 and frac>frac_co and fdr<=fdr_co and t_cat==cat and fo<fo_co:
                                    key = (ss_name, frac_co, 'neg', fdr_co, 'fall_out', fo_co, cat)
                                    counts_dict[key] +=1
                                key = (ss_name, frac_co, 'pos', fdr_co, 'fall_out', fo_co, cat)
                                counts_dict[key] +=0
                                if es>=0 and frac>frac_co and fdr<=fdr_co and t_cat==cat and fo<fo_co:
                                    key = (ss_name, frac_co, 'pos', fdr_co, 'fall_out', fo_co, cat)
                                    counts_dict[key] +=1
            for key in counts_dict.iterkeys():
                ss_name, frac_co, sign, fdr_co, metric, metric_co, cat = key
                count = counts_dict[key]
                lineo = [
                         ss_name,
                         ss_cat,
                         ss_type,
                         frac_co,
                         sign,
                         fdr_co,
                         metric,
                         metric_co,
                         cat,
                         count,
                         ]
                print >>f, '\t'.join(map(str,lineo))
                logging.debug('\t'.join(map(str,lineo)))
    
                
    return 0

if __name__ == '__main__': 
    sys.exit(main())
    
    
    
    
    
    
    

