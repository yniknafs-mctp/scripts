'''
Created on Jul 11, 2014
@author yniknafs
'''

import os
import sys
import argparse
import logging
import collections


'''
takes the SSEA query results for the breast subtypes and nominates a list of lncRNAs associated with each subtype
'''

# pam50_sets = [
#                 'brca_pam_subtype_basal',
#                 'brca_pam_subtype_her2',
#                 'brca_pam_subtype_luma',
#                 'brca_pam_subtype_lumb',
#                 'brca_pam_subtype_normal'
#               ]

pam50_sets = [
                'pam_basal',
                'pam_her2',
                'pam_luma',
                'pam_lumb',
                'pam_normal'
              ]
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("ssea_query")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    meta = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    nomination = collections.defaultdict(lambda: collections.defaultdict(lambda: [0,0]))
    fileh = open(args.ssea_query)
    header = fileh.next().strip().split('\t')
    cn_dict = collections.defaultdict(lambda: [0,0])
    ct_dict = collections.defaultdict(lambda: [0,0])
    #loop through file and store info
    test = []
    for line in fileh: 
        line = line.strip().split('\t')
        t_id = line[header.index('transcript_id')]
        g_id = line[header.index('gene_id')]
        tcat = line[header.index('tcat')]
        set_name = line[header.index('ss_compname')]
        frac = float(line[header.index('frac')])
        fdr = float(line[header.index('fdr')])
        prec = float(line[header.index('prec')])
        fpr = float(line[header.index('fpr')])
        if set_name == 'cancer_versus_normal_breast':
            cn_dict[t_id] = [frac, fdr]
            continue
        if set_name == 'cancer_type_breast_carcinoma':
            ct_dict[t_id] = [frac, fdr]
            continue
        meta[t_id]['g_id'] = g_id
        meta[t_id]['tcat'] = tcat
        nomination[t_id][set_name] = [frac, fdr, prec, fpr]
    
    
    pos_candidates = []
    neg_candidates = []
    herm_candidates = []
    fdr_thresh = 1e-3
    frac_thresh = .75
    prec_thresh = .75
    #check each transcript for significant association with only one subtype 
    for t_id in nomination.iterkeys():
        pos_check = []
        neg_check = []
        cn_frac, cn_fdr = cn_dict[t_id]
        ct_frac, ct_fdr = ct_dict[t_id]
        for xset in nomination[t_id]:
            frac, fdr, prec, fpr = nomination[t_id][xset]
            if frac > 0:
                if prec > prec_thresh:
                    pos_check.append(xset)
            if frac < 0:
                if prec > prec_thresh:
                    neg_check.append(xset)
#             if frac > 0:
#                 if abs(frac) > frac_thresh and fdr < fdr_thresh and prec > prec_thresh:
#                     pos_check.append(xset)
#             if frac < 0:
#                 if abs(frac) > frac_thresh and fdr < fdr_thresh and prec > prec_thresh:
#                     neg_check.append(xset)
            if t_id == 'T036135':
                logging.debug(xset)
                logging.debug('\t'.join(map(str, [frac, fdr, prec, fpr])))
        if len(pos_check) == 1: 
            set_chosen = pos_check[0]
            up_out = [t_id, meta[t_id]['g_id'],meta[t_id]['tcat'], set_chosen, 
                      nomination[t_id][set_chosen][0], nomination[t_id][set_chosen][1],
                      nomination[t_id][set_chosen][2], nomination[t_id][set_chosen][3],  
                      cn_frac, cn_fdr, ct_frac, ct_fdr]
            pos_candidates.append(up_out)
        if len(neg_check) == 1:
            set_chosen = neg_check[0] 
            dn_out = [t_id, meta[t_id]['g_id'],meta[t_id]['tcat'], set_chosen, 
                      nomination[t_id][set_chosen][0], nomination[t_id][set_chosen][1], 
                      nomination[t_id][set_chosen][2], nomination[t_id][set_chosen][3], 
                      cn_frac, cn_fdr, ct_frac, ct_fdr]
            neg_candidates.append(dn_out)
#         if len(pos_check) == 1 and len(neg_check) == 1:
#             set_chosen = pos_check[0]
#             herm_candidates.append([t_id, xset, nomination[t_id][xset][0], nomination[t_id][xset][1]])
#             logging.debug('Positive and Negative hit: %s' %t_id) 
        
    headero = [
               'transcript_id',
               'gene_id',
               'tcat',
               'subtype',
               'frac',
               'fdr',
               'prec',
               'fpr',
               'cn_frac',
               'cn_fdr',
               'ct_frac',
               'ct_fdr'
               ]
    print '\t'.join(headero)
    for cand in pos_candidates: 
        print '\t'.join(map(str, cand))
#     for cand in neg_candidates: 
#         print '\t'.join(map(str, cand))    
    logging.debug(len(pos_candidates))
    return 0

if __name__ == '__main__': 
    sys.exit(main())
