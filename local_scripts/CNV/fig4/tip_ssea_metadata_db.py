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
import shelve

'''
takes a transcript CNV file and obtains transcript
metadata and SSEA stats for the trans
in the peaks
'''
'''
input: transcript file (2col: trans_name, peak_id)
careful with which database you connect to 
'''


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

def db_connect_old(name, host):
    logging.info('connecting to %s database on mongo server: %s' % (name, host))
    client = pymongo.MongoClient(host)
    db = client[name]
    transcripts = db['metadata']
    samples = db['samples']
    sample_sets = db['sample_sets']
    configs = db['configs']
    results = db['reports']
    hists = db['hists']
    merged = db['merged']
    colls = {'transcripts':transcripts, 'samples':samples, 'sample_sets':sample_sets, 
             'configs':configs, 'results':results, 'hists':hists, 'merged':merged}
    return colls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts_file")
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name for ssea run (will be name of database)')
    parser.add_argument("--host", dest = 'host',
                        default = '172.20.54.118:12345',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    #import new metadata
    ####THIS IS THE DATABASE CONTAINING THE UPDATED TRANSCRIPT METADATA#####
    colls_new = db_connect('compendia', '172.20.54.118:27017')
    new_meta = colls_new['transcripts']
    
    #import new metadata
    #####THIS IS THE DATABASE WITH THE SSEA RESULTS###############
    colls = db_connect(args.name, args.host)
    results = colls['results']
    ss = colls['sample_sets']
    transcripts = colls['transcripts']
    logging.info('Parsing transcript file')
    len_file = open(args.transcripts_file,'r')
    tot = len(len_file.readlines())
    i = 0    

    fileo_peak = open('tip_processed_v2.tsv', 'w')
    headers = ['peak_id', 'peak_type', 't_id','gene_id', 'es', 'nes', 'p_value', 'fdr', 'category']
    fileo_peak.write('\t'.join(headers)+'\n')
    tfile = open(args.transcripts_file, 'r')
    for line in tfile:
        i+=1
        sys.stdout.write(("\r completed %d/%d lines of file" % (i, tot)))
        line = line.strip().split('\t')
        t_name = line[0]
        peak_id = line[1]
        type = peak_id.split('-')[1]
        
        #get the sample set id for the peak
        ss_id = ss.find_one({'name':peak_id}, ['_id'])
        if ss_id == None:
            logging.debug("peak %s skipped, not in analysis" % peak_id)
            continue
        ss_id = int(ss_id['_id'])
        
        #get the transcript id 
        t_id = transcripts.find_one({'name':t_name}, ['_id'])
        t_meta = new_meta.find_one({'transcript_id': t_name}, ['transcript_category', 'gene_id'])
        if t_id == None or t_meta == None:
            logging.debug("transcript %s skipped, not in analysis" % t_name)
            continue
        t_id = int(t_id['_id'])
        t_cat = t_meta['transcript_category']
        gene_id = t_meta['gene_id']
        result = results.find_one({'t_id':t_id, 'ss_id':ss_id}, ['es','nes','ss_fdr_q_value','nominal_p_value'])
        nes = float(result['nes'])
        es = float(result['es'])
        fdr = result['ss_fdr_q_value']
        p = result['nominal_p_value']        

        lineo = [peak_id, type, t_name, gene_id, es, nes, p, fdr, t_cat]
        fileo_peak.write('\t'.join(map(str, lineo))+'\n')
    fileo_peak.close()
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
