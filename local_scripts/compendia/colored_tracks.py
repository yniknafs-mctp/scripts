'''
Created on Mar 3, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import json
import collections

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("bed_file")
    parser.add_argument("result_dir")
    parser.add_argument("meta_file")
    parser.add_argument("description")
    parser.add_argument("fdr")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    FDR_THRESH = float(args.fdr)
        
    if args.result_dir.endswith('/'):
        args.result_dir = args.result_dir[:-1] 
    
    #read transcript metadata file to create dictionary for 
    logging.info('reading metadata')
    meta_file = open(args.meta_file)
    id_dict = {}
    id_dict_rev = {}
    for line in meta_file: 
        line_json = json.loads(line.strip())
        id_dict[line_json['transcript_id']] = line_json['_id']
        id_dict_rev[line_json['_id']] = line_json['transcript_id']
    
    
    i=0
    tot=382000
    result_json = os.path.join(args.result_dir, 'results.json')
    es_dict = {} 
    fdr_dict = {}
    frac_dict = {}
    for line in open(result_json): 
        line = line.strip()
        d = json.loads(line)
        t_id = d['t_id']
        i+=1
        if i%5000 == 0:
            logging.debug(("completed %d/%d transcripts" % (i, tot)))
        t_name = id_dict_rev[t_id]
        es = d['es']
        fdr = d['ss_fdr_q_value']
        frac = d['ss_frac']
        es_dict[t_name] = float(es)
        fdr_dict[t_name] = float(fdr)
        frac_dict[t_name] = float(frac)
    
    peak_name = os.path.basename(args.result_dir)
    print "track name=\"" + peak_name + '\" description=\"' + args.description + '\" itemRgb=\"On\"'
    
    t_set = set(es_dict.keys())
    
    result_file = open(args.bed_file)
    for line in result_file: 
        line = line.strip().split('\t')
        #inlcude line below if gene_id is in bed file
#         t_name = line[3].split('|')[1]
        t_name = line[3]
        if t_name not in t_set:
            continue
        fdr = fdr_dict[t_name]
        frac = frac_dict[t_name]
        score = abs(frac)*1000
        
        score_factor = 1-(-1*np.log10(1+1e-3 - abs(frac)))/3.0
        if score_factor <0:
            logging.debug('score_factor: %f, frac: %f' % (score_factor, frac))
        color_score = int(score_factor*255)
        gray = [128,128,128]
        if fdr >= FDR_THRESH: 
            rgb = ','.join(map(str,gray))
        else: 
            if frac >= 0:
                rgb = ','.join(map(str,[255,color_score,color_score]))
            else: 
                rgb = ','.join(map(str,[color_score,color_score,255]))
        line[8]=rgb
        line[4]=score
        print '\t'.join(map(str,line))
        
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
