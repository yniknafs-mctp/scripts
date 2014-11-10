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


gwas_fields = [
               'cancer_snp',
               'overlaps_compendia_exon',
               'Disease/Trait',
               'overlaps_ref_gene',
               'Reported Gene(s)',
               'Tissue'
               ]

GWAS_DIST_THRESH = 25000
FDR_THRESH = 0.01
RANK_THRESH = 380000

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_proc_file")
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
    
    #read GWAS file and save fields to report 
    logging.info('Reading GWAS catalog file')
    gwas_file = open(args.gwas_proc_file)
    gwas_header = gwas_file.next().strip().split("\t")
    gwas_dict = collections.defaultdict(lambda: {})
    for line in gwas_file:
        line = line.strip().split('\t')
        snp_id = line[gwas_header.index('SNPs')]
        for field in gwas_fields: 
            gwas_dict[snp_id][field] = line[gwas_header.index(field)]
    
    # read the sets_type file to make a set_type dict
    logging.info('Reading set type file')
    set_type_file = open(args.set_type_file)
    set_type_header = set_type_file.next().strip().split("\t")
    set_type_dict = {}
    for line in set_type_file:
        line = line.strip().split('\t')
        set_name = line[set_type_header.index('set_name')]
        set_type = line[set_type_header.index('tissue')]
        set_type_dict[set_name] = set_type
    
    colls = db_connect_colls(args.name, args.host)
    transcripts = colls['transcripts']
    sample_sets = colls['sample_sets']
    results = colls['results']

    colls_meta = db_connect_colls(args.name, '172.20.54.118')
    trans_meta = colls_meta['transcripts']

    # query the transcripts database to return transcripts near/on GWAS snps
    gwas_transcripts_query = transcripts.find({'gwas_snp_dist': {'$lte':GWAS_DIST_THRESH}})
    gwas_transcripts = []
    for transcript in gwas_transcripts_query: 
        t_id = transcript['_id']
        snp_id = transcript['gwas_snp_id']
        if gwas_dict[snp_id]!= {}:
            ref_overlap = gwas_dict[snp_id]['overlaps_ref_gene']
            cancer_snp = gwas_dict[snp_id]['cancer_snp']
            reported = gwas_dict[snp_id]['Reported Gene(s)']
            
            if ref_overlap != '1' and cancer_snp=='1': #and (reported == 'intergenic' or reported == 'NR'):
                gwas_transcripts.append(t_id)

    j=0
    tot = len(gwas_transcripts)
    headero = [
               't_name',
               't_cat',
               'ss_name',
               'es',
               'nes',
               'fdr',
               'rank',
               'snp_id',
               'snp_dist',
               'disease',
               'reported_genes',
               'snp_disease_type',
               'set_type'
               ]
    print '\t'.join(headero)
    
    for transcript in gwas_transcripts:
#         hint_index = [('t_id', pymongo.ASCENDING),
#                       ('ss_fdr_q_value', pymongo.ASCENDING)
#                       ]
        results_query = results.find({'ss_rank': {'$lt': RANK_THRESH},
                                      'ss_fdr_q_value': {'$lt': FDR_THRESH},
                                     't_id': transcript}, ['nes', 'ss_fdr_q_value', 'ss_rank', 'es', 'ss_id'])                            
        j+=1
        logging.info('finished %d/%d transcripts' % (j, tot))
        for result in results_query: 
            trans_query = transcripts.find_one({'_id': transcript})
            t_name = trans_query['transcript_id']
            trans_meta_query = trans_meta.find_one({'transcript_id': t_name})
            logging.debug(trans_meta_query)
            if trans_meta_query:
                t_cat = trans_meta_query['transcript_category']
            else: 
                t_cat = 'NULL'
            snp_id = trans_query['gwas_snp_id']
            snp_dist = trans_query['gwas_snp_dist']
            disease = gwas_dict[snp_id]['Disease/Trait']
            reported_gene = gwas_dict[snp_id]['Reported Gene(s)']
            snp_tissue = gwas_dict[snp_id]['Tissue']
            ss_id = result['ss_id']
            ss_name = sample_sets.find_one({'_id': ss_id})['name']
            ss_type = set_type_dict[ss_name]
            nes = result['nes']
            fdr = result['ss_fdr_q_value']
            es = result['es']
            rank = result['ss_rank']
               
            lineo = [t_name, t_cat, ss_name, es, nes, fdr, rank, snp_id, snp_dist, disease, reported_gene, snp_tissue, ss_type]
            print '\t'.join(map(str, lineo))



    
    
    
    logging.info('Running main script')
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
