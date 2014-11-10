'''
Created on Jan 28, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import pymongo
import collections

'''
Takes the GWAS processed file and connects to mongodb 
with all the SSEA results. Finds transcripts near GWAS SNPs
and reports SSEA results for those transcripts
'''


#function to connect to database and return all collections
def db_connect_colls(name, host):
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


GWAS_DIST_THRESH = 25000

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_overlap_file")
    parser.add_argument("set_type_file")
    parser.add_argument("-n", "--name", dest = 'name',
                        default = 'compendia',
                        help = 'name of mongo database')
    parser.add_argument("--host", dest = 'host',
                        default = '172.20.54.118:12345',
                        help = 'name of mongodb server to connect to')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    # read the sets_type file to make a set_type dict
    logging.info('Reading set type file')
    set_type_file = open(args.set_type_file)
    set_type_header = set_type_file.next().strip().split("\t")
    type_sets = collections.defaultdict(lambda: [])
    complex_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: {}))
    for line in set_type_file:
        line = line.replace('\n', '').split('\t')
        set_name = line[set_type_header.index('set_name')]
        set_type = line[set_type_header.index('tissue')]
        set_ca = line[set_type_header.index('ca_code')]
        set_type_code = line[set_type_header.index('type_code')]
        if set_ca and set_type_code: 
            type_sets[set_type].append(set_name)
            complex_dict[set_type][set_ca][set_type_code] = set_name
    

    #read GWAS file and save fields to report 
    logging.info('Reading GWAS file')
    gwas_file = open(args.gwas_overlap_file)
    gwas_header = gwas_file.next().strip().split("\t")
    snp_type_dict = {}
    snp_reported_dict = {}
    for line in gwas_file:
        line = line.strip().split('\t')
        snp_type = line[gwas_header.index('Tissue')]
        snp_id = line[gwas_header.index('SNPs')]
        snp_type_dict[snp_id] = snp_type
        reported = line[gwas_header.index('Reported Gene(s)')]
        snp_reported_dict[snp_id] = reported
        
    colls = db_connect_colls(args.name, args.host)
    transcripts = colls['transcripts']
    sample_sets = colls['sample_sets']
    results = colls['results']

    colls_meta = db_connect_colls(args.name, '172.20.54.118')
    trans_meta = colls_meta['transcripts']

    #print header
    headero = [
               't_name',
               't_cat',
               'snp_id',
               'snp_dist',
               'cancer_type',
               'reported_genes',
               'ca_nes',
               'ca_fdr',
               'ca_es',
               'ca_rank',
               'type_nes',
               'type_fdr',
               'type_es',
               'type_rank'
               ]
    print '\t'.join(headero)

    # query the transcripts database to return transcripts near/on GWAS snps
    gwas_transcripts_query = transcripts.find({'gwas_snp_dist': {'$lte':GWAS_DIST_THRESH}})
    gwas_transcripts = []
    cancer_snps = set(snp_type_dict.keys())
    tot = gwas_transcripts_query.count()
    gwas_transcripts_query = transcripts.find({'gwas_snp_dist': {'$lte':GWAS_DIST_THRESH}})
    j = 0
    for transcript in gwas_transcripts_query: 
        j+=1
        if j%2500==0:
            logging.debug('finished %d/%d transcripts' % (j, tot))
        t_id = transcript['_id']
        snp_id = transcript['gwas_snp_id']
        if snp_id not in cancer_snps: 
            continue
        else: 
            reported = snp_reported_dict[snp_id]
            snp_type = snp_type_dict[snp_id]
            t_name = transcript['transcript_id']
            trans_meta_query = trans_meta.find_one({'transcript_id': t_name})
            t_cat = trans_meta_query['transcript_category']
            snp_dist = transcript['gwas_snp_dist']
            for cancer in complex_dict[snp_type].iterkeys():
                cancer_dict = complex_dict[snp_type][cancer]
                ca_set_name = cancer_dict['cancer']
                type_set_name = cancer_dict['type']
                ca_set_id = sample_sets.find_one({'name': ca_set_name})['_id']
                type_set_id = sample_sets.find_one({'name': type_set_name})['_id']
                #get results for the ca set
                ca_results = results.find_one({'t_id': t_id, 'ss_id': ca_set_id})
                ca_nes = ca_results['nes']
                ca_fdr = ca_results['ss_fdr_q_value']
                ca_es = ca_results['es']
                ca_rank = ca_results['ss_rank']
                #now results for the type set
                type_results = results.find_one({'t_id': t_id, 'ss_id': type_set_id})
                type_nes = type_results['nes']
                type_fdr = type_results['ss_fdr_q_value']
                type_es = type_results['es']
                type_rank = type_results['ss_rank']
                
                #print line
                lineo = [t_name, t_cat, snp_id, snp_dist, cancer, reported, 
                         ca_nes, ca_fdr, ca_es, ca_rank,
                         type_nes, type_fdr, type_es, type_rank]
                print '\t'.join(map(str, lineo))

    
    
    
    logging.info('Running main script')
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
