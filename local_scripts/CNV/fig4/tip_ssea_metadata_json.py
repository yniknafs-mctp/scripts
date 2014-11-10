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

'''
takes a transcript CNV file and obtains transcript
metadata and SSEA stats for the trans
in the peaks
'''
'''
input: transcript file (2col: trans_name, peak_id)
will loop through all the json files 
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
    parser.add_argument("dels_dir")
    parser.add_argument("amps_dir")
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name for ssea run (will be name of database)')
    parser.add_argument("--host", dest = 'host',
                        default = '172.20.54.118',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #make a dict of file location for the results json for differnet peaks
    peak_file_dict = {}
    #amps first
    if args.amps_dir != 'null':
        amp_file_path = os.path.join(args.amps_dir, 'sample_sets.tsv')
        amp_file = open(amp_file_path)
        header = amp_file.next().strip().split('\t')
        for line in amp_file:
            line = line.strip().split('\t')
            dirname = line[header.index('dirname')]
            result_json = os.path.join(args.amps_dir, dirname, 'results.json')
            peak_name = line[header.index('name')]
            peak_type = peak_name.split('-')[1]
            if peak_type == 'Amp': 
                peak_file_dict[peak_name] = result_json
    else: 
        logging.debug("Not parsing amplification peaks")
    #then dels
    if args.dels_dir != 'null': 
        del_file_path = os.path.join(args.dels_dir, 'sample_sets.tsv')
        del_file = open(del_file_path)
        header = del_file.next().strip().split('\t')
        for line in del_file:
            line = line.strip().split('\t')
            dirname = line[header.index('dirname')]
            result_json = os.path.join(args.dels_dir, dirname, 'results.json')
            peak_name = line[header.index('name')]
            peak_type = peak_name.split('-')[1]
            if peak_type == 'Del': 
                peak_file_dict[peak_name] = result_json
    else: 
        logging.debug("Not parsing deletion peaks")
    
    fields_trans = ['chrom',
                'start',
                'end',
                'transcript_category',
                'transcript_id',
                'gene_id',
                '_id']
    colls = db_connect(args.name, args.host)
    transcripts = colls['transcripts']
    #parse through transcript metadata and build dict to be used during merge
    trans_dict = collections.defaultdict(lambda: {})
    trans_dict_by_id = collections.defaultdict(lambda: {})
    logging.info('Parsing through transcript metadata')
    tot = transcripts.find().count()
    i = 0
    for x in transcripts.find():
        i+=1
        if (i % 50000) == 0:
            logging.debug('Finished %d/%d' % (i, tot))
            
        key = x['transcript_id']
        key2 = x['_id']
        #create a dict placeholder for this _id element
        id_dict = {}
        for field in fields_trans:
            id_dict[field] = x[field]
        #create a combined locus and strand field
        trans_dict[key] = id_dict
        trans_dict_by_id[key2] = id_dict
    
    
    fileo_peak = open('tip_processed_v3.tsv', 'w')
    headers = ['peak_id', 'peak_type', 't_id','gene_id', 'es', 'nes', 'p_value', 'fdr', 'category']
    fileo_peak.write('\t'.join(headers)+'\n')
    tip_dict = collections.defaultdict(lambda: [])
    tfile = open(args.transcripts_file, 'r')
    tot = 0
    i = 0
    for line in tfile:
        line = line.strip().split('\t')
        t_name = line[0]
        peak_id = line[1]
        if trans_dict[t_name] == {}:
            continue
        t_id = trans_dict[t_name]['_id']
        type = peak_id.split('-')[1]
        tip_dict[peak_id].append(t_id)
        if args.amps_dir != 'null':
            if type == 'Amp':
                tot+=1
        if args.dels_dir != 'null':
            if type == 'Del':
                tot+=1
    
    for peak in peak_file_dict.iterkeys():
        print peak
    
    for peak in tip_dict.iterkeys():
        #open the sample set file of this peak
#         peak = peak.replace('-','_')
        type = peak.split('-')[1]
        if args.amps_dir == 'null':
            if type == 'Amp':
                continue
        if args.dels_dir == 'null':
            if type == 'Del':
                continue
        
        ss_file_path = peak_file_dict[peak]
        ss_file = open(ss_file_path) 
        tips = set(tip_dict[peak])
        type = peak.split('-')[1]
        
        for line in ss_file: 
            line = line.strip()
            d = json.loads(line)
            t_id = d['t_id']
            if t_id in tips:
                i+=1
                sys.stdout.write(("\r completed %d/%d transcripts" % (i, tot)))
                t_meta = trans_dict_by_id[t_id]
                t_name = t_meta['transcript_id']
                gene_id = t_meta['gene_id']
                chrom = t_meta['chrom']
                start = t_meta['start']
                end = t_meta['end']
                cat = t_meta['transcript_category']
                nes = d['nes']
                fdr = d['ss_fdr_q_value']
                p = d['nominal_p_value']
                es = d['es']
                lineo = [peak, type, t_name, gene_id, es, nes, p, fdr, cat]
                fileo_peak.write('\t'.join(map(str, lineo))+'\n')
    fileo_peak.close()
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
