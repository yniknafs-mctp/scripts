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



'''
reads a tip_file and reports ssea info for each tip
'''


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("tip_file")
    parser.add_argument("result_dir")
    parser.add_argument("meta_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read transcript metadata file to create dictionary for metadata
    logging.info('reading transcript metadata')
    meta_file = open(args.meta_file)
    meta_dict = {}
    id_dict = {}
    id_dict_rev = {}
    gene_id_dict = {}
    for line in meta_file: 
        line_json = json.loads(line.strip())
        meta_dict[line_json['_id']] = line_json['tcat']
        id_dict[line_json['transcript_id']] = line_json['_id']
        id_dict_rev[line_json['_id']] = line_json['transcript_id']
        gene_id_dict[line_json['_id']] = line_json['gene_id']
    
    
    # read tip file and make a dict. key = peak, values = tips
    logging.info('reading tip_file')
    tfile = open(args.tip_file)
    tot = 0
    tip_dict = collections.defaultdict(lambda: [])
    tip_set = set(id_dict.keys())
    for line in tfile:
        tot+=1
        line = line.strip().split('\t')
        t_name = line[0]
        peak_id = line[1]
        if t_name not in tip_set:
            continue
        t_id = id_dict[t_name]
        tip_dict[peak_id].append(t_id)

    #loop through all the peaks in the peak dict, open results file, and get data    
    i=0
    logging.info('parsing peaks')
    for peak in tip_dict.iterkeys():
        logging.debug('Processing peak: %s', peak)
        #open the sample set file of this peak
        type = peak.split('-')[1]
        peak_dir = peak.replace('-','_').replace('.','_').lower()
        result_json = os.path.join(args.result_dir, peak_dir, 'results.json') 
        tips = set(tip_dict[peak])
        for line in open(result_json): 
            line = line.strip()
            d = json.loads(line)
            t_id = d['t_id']
            if t_id in tips:
                i+=1
                if i%50 == 0:
                    logging.debug(("completed %d/%d transcripts" % (i, tot)))
                t_name = id_dict_rev[t_id]
                cat = meta_dict[t_id]
                g_id = gene_id_dict[t_id]
                nes = d['nes']
                fdr = d['ss_fdr_q_value']
                p = d['nominal_p_value']
                es = d['es']
                lineo = [peak_dir, type, t_name, g_id, es, nes, fdr, cat]
                print '\t'.join(map(str,lineo))
                
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
