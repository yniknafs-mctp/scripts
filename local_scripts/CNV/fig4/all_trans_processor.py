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
Parses through JSON files in an SSEA folder
will append transcript metadata to each result file 
and reprint it. 

This is used for attaing the NES scores for all transcripts
everywhere for each peak
'''
'''
input: transcript CNV file (2 cols: trans_name, peak_name)
        directory to find all sample sets for that CNV analysis 
'''


fields_trans = ['chrom',
                'start',
                'end',
                'transcript_category',
                'transcript_id',
                'gene_id']

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
    parser.add_argument("transcript_cnv_file")
    parser.add_argument("sample_sets_dir")
    parser.add_argument("output_dir")
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name for ssea run (will be name of database)')
    parser.add_argument("--host", dest = 'host',
                        default = '172.20.54.118',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Parsing transcript cnv file')
    tfile = open(args.transcript_cnv_file, 'r')
    peaks_set = set()
    for line in tfile:
        line = line.strip().split('\t')
        t_name = line[0]
#         peak_id = line[1].replace('-','_')
#         type = peak_id.split('_')[1]
        peaks_set.add(t_name)
    
    colls = db_connect(args.name, args.host)
    transcripts = colls['transcripts']
    
    #parse through transcript metadata and build dict to be used during merge
    trans_dict = collections.defaultdict(lambda: {})
    logging.info('Parsing through transcript metadata')
    tot = transcripts.find().count()
    i = 0
    for x in transcripts.find():
        i+=1
        if (i % 50000) == 0:
            logging.debug('Finished %d/%d' % (i, tot))
            
        key = x['_id']
        #create a dict placeholder for this _id element
        id_dict = {}
        for field in fields_trans:
            id_dict[field] = x[field]
        #create a combined locus and strand field
        trans_dict[key] = id_dict

    ss_file_dir = os.path.join(args.sample_sets_dir, 'sample_sets.tsv')
    ss_tot = len(open(ss_file_dir).readlines())
    ss_file = open(ss_file_dir)
    ss_header = ss_file.next().strip().split('\t')
    logging.info('Generating tsv files from results jsons')
    j = 0
    for line in ss_file:
        j+=1
        line = line.strip().split('\t')
        dirname = line[ss_header.index('dirname')]
        result_json = os.path.join(args.sample_sets_dir, dirname, 'results.json')
        peak_name = line[ss_header.index('name')]
        peak_type = peak_name.split('-')[1]
        logging.info("Processing peak %d/%d: %s" % (j, ss_tot, peak_name))
        
        tot = len(open(result_json).readlines())
        results = open(result_json)
        
        fileo_name = peak_name + '_peak_member.tsv'
        fileo_dir = os.path.join(args.output_dir, fileo_name)
        fileo = open(fileo_dir, 'w')
        header = ['t_id', 'gene_id', 'chrom', 'start', 'end', 'category', 'peak_name', 'peak_membership', 'peak_type', 'nes', 'fdr', 'p']
        fileo.write('\t'.join(header)+'\n')
        i = 0
        for line in results:
            i+=1
            if (i%10000) == 0: 
                logging.debug("Finished %d/%d transcripts for %s: peak %d/%d" % (i, tot, peak_name, j, ss_tot)) 
            line = line.strip()
            d = json.loads(line)
            t_id = d['t_id']
            t_meta = trans_dict[t_id]
            t_name = t_meta['transcript_id']
            gene_id = t_meta['gene_id']
            chrom = t_meta['chrom']
            start = t_meta['start']
            end = t_meta['end']
            cat = t_meta['transcript_category']
            if t_name in peaks_set:
                in_peak = 1
            else: 
                in_peak = 0
            nes = d['nes']
            fdr = d['ss_fdr_q_value']
            p = d['nominal_p_value']
            
            lineo = [t_name, gene_id, chrom, start, end, cat, peak_name, in_peak, peak_type, nes, fdr, p]
            fileo.write("\t".join(map(str, lineo))+'\n')
        
#     headers = ['peak_id', 't_id', 'nes', 'category', 'peak_type']  
#     fileo.write('\t'.join(headers)+'\n')
#     for key in best_trans.iterkeys(): 
#         tup = best_trans[key]
#         lineo = [key] + tup
#         fileo.write('\t'.join(map(str, lineo))+ '\n')        
#     fileo.close()
    return 0

if __name__ == '__main__': 
    sys.exit(main())
