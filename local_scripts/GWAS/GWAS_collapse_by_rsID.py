import argparse
import logging
import os
import sys
import csv
import time
import collections
import bisect
import operator
import subprocess

'''
Takes the tab delimited  GWAS catalog and collapses rows by SNPid
Some SNPs are in the catalog multiple times. 

Input: tab delimited GWAS catalog
'''

gwas_fields_combine = ['PUBMEDID',
                       'Disease/Trait']
gwas_fields = ['SNPs',
               'Snp_id_current',
               'Chr_id',
               'Chr_pos',
               'Reported Gene(s)',
               'Mapped_gene',
               'Upstream_gene_id',
               'Downstream_gene_id',
               'Upstream_gene_distance',
               'Downstream_gene_distance',
               'Snp_gene_ids',
               'Context',
               'Intergenic']

def main():

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("GWAS_file")
    args = parser.parse_args()
    logging.info("Starting")
    rsID_dict = collections.defaultdict(lambda: {})
    rsID_pubmed_dict = collections.defaultdict(lambda: [])
    rsID_disease_dict = collections.defaultdict(lambda: [])
    with open(args.GWAS_file, 'r') as f:
        headers = f.next().strip().split('\t')
        for row in f:
            row = row.strip().split('\t')
            if row != ['']:
                #sure that for each SNP a single SNPid is listed and a genomic position is listed
                ID = row[headers.index('Snp_id_current')]
                if not ID: 
                    continue
                chrom = row[headers.index('Chr_id')]
                if not chrom: 
                    continue
                #combine PUBMEDID and disease for same SNPS in different studies
                pubmedID = row[headers.index('PUBMEDID')]
                rsID_pubmed_dict[ID].append(pubmedID)
                disease = row[headers.index('Disease/Trait')]
                rsID_disease_dict[ID].append(disease)
                #save all fields for SNP to be printed 
                if rsID_dict[ID] == {}:
                    for field in gwas_fields: 
                        row_field = row[headers.index(field)]
                        rsID_dict[ID][field] = row_field
    #print GWAS results unique to each SNP
    fileo = open('gwas_rsID.tsv','w')
    headero = gwas_fields + gwas_fields_combine
    print >> fileo, '\t'.join(map(str,headero))
    for key in rsID_dict.iterkeys():
        lineo = []
        rsID_vals = rsID_dict[key]
        for field in gwas_fields: 
            field_value = rsID_vals[field]
            lineo.append(field_value)
        pubmedIDs = rsID_pubmed_dict[key]
        lineo.append(','.join(map(str, pubmedIDs)))
        diseases = rsID_disease_dict[key]
        lineo.append(','.join(map(str, diseases)))
        print >> fileo, '\t'.join(map(str, lineo))
    fileo.close()
    
    
            
    return 0



if __name__ == '__main__':

    sys.exit(main())