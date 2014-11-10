'''
Created on Jan 4, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import pymongo
import collections
import json

fields_trans = ['chrom',
                'start',
                'end',
                'transcript_category',
                'transcript_id']

def db_connect(name, host):
    logging.info('connecting to %s database on mongo server: %s' % (name, host))
    client = pymongo.MongoClient(host)
    db = client[name]
    transcripts = db['transcripts']
    samples = db['samples']
    sample_sets = db['sample_sets']
    configs = db['configs']
    results = db['results']
    hists = db['hists']
    merged = db['merged']
    colls = {'transcripts':transcripts, 'samples':samples, 'sample_sets':sample_sets, 
             'configs':configs, 'results':results, 'hists':hists, 'merged':merged}
    return colls

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("conservation_file")
    parser.add_argument("sample_set_file")
    parser.add_argument("output_dir")
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name for ssea run (will be name of database)')
    parser.add_argument("--host", dest = 'host',
                        default = 'localhost:27017',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    

    
    #parse through transcript metadata and build dict to be used during merge
    trans_dict = collections.defaultdict(lambda: 'NA')
    logging.info('Parsing through transcript metadata')
    colls = db_connect(args.name, args.host)
    transcripts = colls['transcripts']
    tot = transcripts.find().count()
    i = 0
    for x in transcripts.find():
        i+=1
        if (i % 50000) == 0:
            logging.debug('Finished %d/%d' % (i, tot))
        value = x['_id']
        key = x['transcript_id']
        trans_dict[key] = value

    logging.info('Parsing through results json for SSEA values')
    tot = len(open(args.sample_set_file).readlines())
    results = open(args.sample_set_file)
    i = 0
    result_dict = collections.defaultdict(lambda: [])
    for line in results:
        i+=1
        if (i%10000) == 0: 
            logging.debug("Finished %d/%d transcripts" % (i, tot)) 
        line = line.strip()
        d = json.loads(line)
        t_id = d['t_id']
        es = d['es']
        nes = d['nes']
        fdr = d['ss_fdr_q_value']
        p = d['nominal_p_value']
        val = [es, nes, fdr, p]
        result_dict[t_id] = val

    fileo_name = 'assembly.metadata.allcan_cons.tsv'
    fileo_dir = os.path.join(args.output_dir, fileo_name)
    fileo = open(fileo_dir, 'w')
    tot = len(open(args.conservation_file).readlines())
    logging.info('Parsing conservation file')
    tfile = open(args.conservation_file, 'r')
    header = tfile.next().strip().split('\t')
    headero = header + ['allcan_es', 'allcan_nes', 'allcan_fdr', 'allcan_p']
    i=0
    print >>fileo, '\t'.join(headero)
    for line in tfile:
        i+=1
        if (i%10000) == 0: 
            logging.debug("Finished %d/%d transcripts" % (i, tot))
        line = line.strip().split('\t')
        t_name = line[0]
        json_t_id = trans_dict[t_name]
        if json_t_id == 'NA':
            ssea_vals = ['NA']*4
            lineo = line + ssea_vals
            print >>fileo, '\t'.join(map(str,lineo))
            continue
        ssea_vals = result_dict[json_t_id]
        lineo = line + ssea_vals
        print >>fileo, '\t'.join(map(str,lineo))
        
    logging.info('Finished')
    return 0

if __name__ == '__main__': 
    sys.exit(main())
