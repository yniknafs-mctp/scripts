'''
Created on Mar 14, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("metadata")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    
    #build metadata dict
    logging.info('reading metadata')
    meta_file = open(args.metadata)
    meta_header = meta_file.next().strip().split('\t')
    meta_dict = {}
    g_dict = {}
    chrom_dict = {}
    start_dict = {}
    end_dict = {}
    status_dict = {}
    genic_dict = {}
    k=0
    for line in meta_file:
        k+=1
        if k%10000==0:
            logging.debug('Parsed %d metadata results' %k) 
        line = line.strip().split('\t')
        t_id = line[meta_header.index("transcript_id")]
        g_id = line[meta_header.index("gene_id")]
        tcat = line[meta_header.index("tcat")]
        chrom = line[meta_header.index("chrom")]
        start = line[meta_header.index("start")]
        end = line[meta_header.index("end")]
        tstatus = line[meta_header.index("tstatus")]
        tgenic = line[meta_header.index("tgenic")]
        meta_dict[t_id] = tcat
        g_dict[t_id] = g_id
        chrom_dict[t_id] = chrom
        start_dict[t_id] = start 
        end_dict[t_id] = end
        status_dict[t_id] = tstatus
        genic_dict[t_id] = tgenic
    
    results_f = open(args.results)
    results_h = results_f.next().strip().split('\t')
    results_h.append('tcat')
    set_dict = collections.defaultdict(lambda: set())
    t_dict = collections.defaultdict(lambda: {})
    for line in results_f:         
        line = line.strip().split('\t')
        set_name = line[results_h.index('ss_compname')]
        t_id = line[results_h.index('transcript_id')]
        es = line[results_h.index('es')]
        fdr = line[results_h.index('fdr')]
        frac = line[results_h.index('frac')]
        t_dict[t_id][set_name] = [es, fdr, frac]
        set_dict[set_name].add(t_id)
    
    vs_norm = set_dict['brca_tnbc_versus_normal']
    vs_non = set_dict['brca_tnbc_vs_nontnbc']
    catype = set_dict['cancer_type_breast_carcinoma']
    union = vs_norm & vs_non & catype
    union_wo_type = vs_norm & vs_non 
    header = [
              'transcript_id',
              'gene_id',
              'transcript_cat',
              'tstatus',
              'tgenic',
              'chrom',
              'start',
              'end',
              'vs_norm_es',
              'vs_norm_fdr',
              'vs_norm_frac',
              'vs_non_es',
              'vs_non_fdr',
              'vs_non_frac',
              'type_es',
              'type_fdr',
              'type_frac',
              'mean_frac_w_type',
              'mean_frac_wo_type'
              ]
    print '\t'.join(header)
    for key in t_dict.iterkeys():
        if key in union: 
            t_id = key
            g_id = g_dict[t_id]
            t_cat = meta_dict[t_id]
            chrom = chrom_dict[t_id]
            start = start_dict[t_id]
            end = end_dict[t_id]
            tstatus = status_dict[t_id]
            tgenic = genic_dict[t_id]
            vs_norm_es, vs_norm_fdr, vs_norm_frac = t_dict[key]['brca_tnbc_versus_normal']
            vs_non_es, vs_non_fdr, vs_non_frac = t_dict[key]['brca_tnbc_vs_nontnbc']
            type_es, type_fdr, type_frac = t_dict[key]['cancer_type_breast_carcinoma']
            avg_w_type = (float(vs_norm_frac) + float(vs_non_frac) + float(type_frac))/3.0
            avg_wo_type = (float(vs_norm_frac) + float(vs_non_frac))/2.0
            lineo = [t_id, g_id, t_cat, 
                     tstatus, tgenic,
                     chrom, start, end,
                     vs_norm_es, vs_norm_fdr, vs_norm_frac, 
                     vs_non_es, vs_non_fdr, vs_non_frac, 
                     type_es, type_fdr, type_frac, 
                     avg_w_type, avg_wo_type]
            print '\t'.join(map(str, lineo))
        elif key in union_wo_type: 
            t_id = key
            g_id = g_dict[t_id]
            t_cat = meta_dict[t_id]
            chrom = chrom_dict[t_id]
            start = start_dict[t_id]
            end = end_dict[t_id]
            tstatus = status_dict[t_id]
            tgenic = genic_dict[t_id]
            vs_norm_es, vs_norm_fdr, vs_norm_frac = t_dict[key]['brca_tnbc_versus_normal']
            vs_non_es, vs_non_fdr, vs_non_frac = t_dict[key]['brca_tnbc_vs_nontnbc']
            avg_w_type = 'NA'
            avg_wo_type = (float(vs_norm_frac) + float(vs_non_frac))/2.0
            lineo = [t_id, g_id, t_cat, 
                     tstatus, tgenic,
                     chrom, start, end,
                     vs_norm_es, vs_norm_fdr, vs_norm_frac, 
                     vs_non_es, vs_non_fdr, vs_non_frac, 
                     'NA', 'NA', 'NA', 
                     avg_w_type, avg_wo_type]
            print '\t'.join(map(str, lineo))
            
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
