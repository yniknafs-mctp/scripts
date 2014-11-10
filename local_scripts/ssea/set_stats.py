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

'''
Takes file containing the path for all set results 
to be used in the set stats report 
'''


header = [
                 'ss_name',
                 'fdr_cutoff',
                 'num_tot',
                 'num_lnc',
                 'num_tucp',
                 'pos_es_max',
                 'pos_es_75', 
                 'pos_es_med',
                 'pos_es_25', 
                 'pos_es_min',
                 'neg_es_max',
                 'neg_es_75', 
                 'neg_es_med',
                 'neg_es_25', 
                 'neg_es_min',
                 'pos_nes_max',
                 'pos_nes_75',
                 'pos_nes_med',
                 'pos_nes_25',
                 'pos_nes_min',
                 'neg_nes_max', 
                 'neg_nes_75', 
                 'neg_nes_med',
                 'neg_nes_25', 
                 'neg_nes_min',

            ]

fdr_cutoffs = [1e-2, 1e-3, 1e-4, 1e-5, 0]

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata_json")
    parser.add_argument("sets_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    #read transcript metadata file to create dictionary for 
    logging.info('reading metadata')
    meta_file = open(args.metadata_json)
    meta_dict = {}
    for line in meta_file: 
        line_json = json.loads(line.strip())
        meta_dict[line_json['_id']] = line_json
    
    print >> sys.stdout, '\t'.join(header)
    # read through SS file to process each SS
    sets_file = open(args.sets_file)
    for ss in sets_file:
        ss = ss.strip()
        #read sample_set json to obtain sample set name
        sample_set = json.load(open(os.path.join(ss, 'sample_set.json')))
        ss_name = sample_set['name']
        logging.info('Processing sample set: %s' % ss_name)
        #reinitialize stat lists
        fdrs = []
        ess = []
        ness = []
        cats = []
        #read through result json to obtain result stats
        results = open(os.path.join(ss, 'results.json'))
        results_read = results.readlines()
        i = 0
        for x in results_read:
            i+=1
            if (i%25000) == 0: 
                logging.debug('Finished %d/%d transcripts' % (i, len(results_read)))
            result = json.loads(x.strip())
            t_id = result['t_id']
            t_cat = meta_dict[t_id]['tcat']
            fdrs.append(result['ss_fdr_q_value'])
            ess.append(result['es'])
            ness.append(result['nes'])
            cats.append(t_cat)
            
        fdrs = np.array(fdrs)
        ess = np.array(ess)
        ness = np.array(ness)
        cats = np.array(cats)
        
        neg_es = ess[ess<=0]
        pos_es = ess[ess>0]
        neg_nes = ness[ess<=0]
        pos_nes = ness[ess>0]
        
        for cutoff in fdr_cutoffs: 
            lt = fdrs<=cutoff            
            num_lt = sum(lt)
            lnc_lt = sum(cats[lt]=='lncrna')
            tucp_lt = sum(cats[lt]=='tucp')
            
            pos_es_max = pos_es.max()
            pos_es_min = pos_es.min()
            pos_es_25 = np.percentile(pos_es, 25)
            pos_es_75 = np.percentile(pos_es, 75)
            pos_es_med = np.median(pos_es)
             
            neg_es_max = neg_es.max()
            neg_es_min = neg_es.min()
            neg_es_25 = np.percentile(neg_es, 25)
            neg_es_75 = np.percentile(neg_es, 75)
            neg_es_med = np.median(neg_es)
             
            pos_nes_max = pos_nes.max()
            pos_nes_min = pos_nes.min()
            pos_nes_25 = np.percentile(pos_nes, 25)
            pos_nes_75 = np.percentile(pos_nes, 75)
            pos_nes_med = np.median(pos_nes)
             
            neg_nes_max = neg_nes.max()
            neg_nes_min = neg_nes.min()
            neg_nes_25 = np.percentile(neg_nes, 25)
            neg_nes_75 = np.percentile(neg_nes, 75)
            neg_nes_med = np.median(neg_nes)
            
            lineo = [
                     ss_name,
                     cutoff,
                     num_lt,
                     lnc_lt,
                     tucp_lt,
                     pos_es_max,
                     pos_es_75, 
                     pos_es_med,
                     pos_es_25, 
                     pos_es_min,
                     neg_es_max,
                     neg_es_75, 
                     neg_es_med,
                     neg_es_25, 
                     neg_es_min,
                     pos_nes_max,
                     pos_nes_75,
                     pos_nes_med,
                     pos_nes_25,
                     pos_nes_min,
                     neg_nes_max, 
                     neg_nes_75, 
                     neg_nes_med,
                     neg_nes_25, 
                     neg_nes_min
                     ]
        
            print >> sys.stdout, '\t'.join(map(str,lineo))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
    
    
    
    
    
    
    

