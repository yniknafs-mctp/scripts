'''
Created on May 23, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections


'''
take gwas snp catalog, and the file with TIDs matched to any significant associations,
and the file that keys in the sample set name to a GWAS flagged tissue type 
and report any TID that has a significant association 
'''


def compnamer(compname):
    if compname.startswith('cancer_type'):
        pre = 'Lineage Specificty'
        mid = ' '.join(compname.split('cancer_type')[1].split('_'))
        out = pre + mid.upper()
    if compname.startswith('normal_type'):
        pre = 'Lineage Specificty'
        mid = ' '.join(compname.split('normal_type')[1].split('_'))
        out = pre + mid.upper()
    if compname.startswith('cancer_versus_normal'):
        pre = 'Cancer vs Normal'
        mid = ' '.join(compname.split('cancer_versus_normal')[1].split('_'))
        out = pre + mid.upper()
    return out


def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_id")
    parser.add_argument("gwas_catalog")
    parser.add_argument("ssea_results")
    parser.add_argument("gwas_overlap")
    parser.add_argument("meta")
    parser.add_argument('--frac', dest='frac', type=float, default=0.90)
    parser.add_argument('--dist', dest='dist', type=int, default=0)
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read transcript metadata
    logging.debug("reading transcript metadata")
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_dict = collections.defaultdict(lambda: 'null')
    for line in meta_fh:
        lineo = line.strip().split('\t')
        t_id = lineo[meta_header.index('transcript_id')]        
        meta_dict[t_id] = lineo
    
    #read transcript metadata to key each snp into the neareast compendia transcripts
    logging.debug('Keying each snp to nearest transcript')
    meta_fh = open(args.gwas_overlap)
    SNP_COL = 3
    TID_COL = 7
    DIST_COL = 16
    snp_dict = collections.defaultdict(lambda: [])
    snp_dist_dict = collections.defaultdict(lambda: [])
    for line in meta_fh: 
        line = line.strip().split('\t')
        t_id = line[TID_COL]
        snp_id = 'rs'+line[SNP_COL]
        snp_dist = int(line[DIST_COL])
        if snp_dist <= args.dist:
            snp_dict[t_id].append(snp_id)
            snp_dist_dict[t_id].append(snp_dist)
    
    #read the category file to make a dictionary to translate ssea directories into GWAS tissue names
    gwas_id_fh = open(args.gwas_id)
    gwas_ssea_dict = collections.defaultdict(lambda: 'null')
    for line in gwas_id_fh:
        gwas, ssea = line.strip().split('\t')
        gwas_ssea_dict[ssea] = gwas
      
    # read the snp catalog (for all the intergenic cancer snps) file to key each snp into a tissue type (for cancer)
    logging.debug('Reading SNP catalog to assign tissue to each snp')
    snp_fh = open(args.gwas_catalog)
    snp_header = snp_fh.next().strip().split('\t')
    snp_tissue_dict = collections.defaultdict(lambda: 'null')
    snp_disease_dict = collections.defaultdict(lambda: 'null')
    snp_pubmed_dict = collections.defaultdict(lambda: 'null')
    for line in snp_fh: 
        line = line.strip().split('\t')
        snp_id = line[snp_header.index("SNPs")]
        tissue = line[snp_header.index("Tissue")]
        snp_tissue_dict[snp_id] = tissue
        disease = line[snp_header.index("Disease/Trait")]
        snp_disease_dict[snp_id] = disease
        pubmed = line[snp_header.index("PUBMEDID")]
        snp_pubmed_dict[snp_id] = pubmed
    
    #print header for the table
    headero = [
               'snp_id',
               'snp_disease',
               'snp_pubmed',
               't_id',
               'g_id',
               'func_name',
               'loc',
               'tcat',
               'snp_dist',
               'set_name',
               'frac'
               ]
    print '\t'.join(headero)
    
    #read the ssea result file to key each transcript with a significant gwas association 
    logging.debug("Parsing SSEA results to find matches")
    ssea_fh = open(args.ssea_results)
    ssea_header = ssea_fh.next().strip().split("\t") 
    for line in ssea_fh: 
        line = line.strip().split('\t')
        t_id = line[ssea_header.index('transcript_id')]
        setname = line[ssea_header.index('ss_compname')]
        frac = abs(float(line[ssea_header.index('frac')]))
        frac_noabs = float(line[ssea_header.index('frac')])
        if frac <= args.frac:
            continue
        gwas_name = gwas_ssea_dict[setname]
        snps = snp_dict[t_id]
        snp_dists = snp_dist_dict[t_id]
        for x in xrange(len(snps)):
            snp = snps[x]
            snp_dist = snp_dists[x]
            snp_tissue = snp_tissue_dict[snp]
            snp_disease = snp_disease_dict[snp]
            snp_pubmed = snp_pubmed_dict[snp]
            if snp_tissue != 'null' and snp_tissue !='NA':
                if snp_tissue == gwas_name:
                    if meta_dict[t_id] == 'null':
                        continue
                    g_id = meta_dict[t_id][meta_header.index('gene_id')]
                    chrom = meta_dict[t_id][meta_header.index('chrom')]
                    start = meta_dict[t_id][meta_header.index('start')]
                    end = meta_dict[t_id][meta_header.index('end')]
                    strand = meta_dict[t_id][meta_header.index('strand')]
                    tcat = meta_dict[t_id][meta_header.index('tcat')]
                    func_name = meta_dict[t_id][meta_header.index('func_name')]
                    if tcat not in ['lncrna', 'tucp']:
                        continue
                    loc = chrom + ':' + start + '-' + end + '(' + strand + ')'
                    lineo = [snp, snp_disease, snp_pubmed, t_id, g_id, func_name, loc, tcat, snp_dist, compnamer(setname), frac_noabs]
                    print '\t'.join(map(str,lineo))
        
        
    
        
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
