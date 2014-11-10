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
Takes the gwas_rsID.tsv file and other bedtools output files 
adds info to the gwas_rsID file

input: unique SNP tsv (gwas_rsID.tsv from GWAS_collapse_by_rsID.py)
'''

def main():

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("gwas_rsID_file")
    parser.add_argument("bedtools_trans")
    parser.add_argument("bedtools_exons")
    parser.add_argument("bedtools_trans_ref")
    parser.add_argument("bedtools_exons_ref")
    parser.add_argument("bedtools_intergenic")
    parser.add_argument("cancer_catalog")
    args = parser.parse_args()
    logging.info("Starting")
    
    #open bedtools overlap for whole transcript
    #store value for whether or not the SNP overlaps a compendia gene    
    RSIDCOL = 3
    in_compendia = collections.defaultdict(lambda: 0)
    with open(args.bedtools_trans, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            rdID = line[RSIDCOL]
            in_compendia[rdID] = 1
            
    #open bedtools overlap for exon overlap
    #store value for whether or not the SNP overlaps a compendia exon    
    in_exon = collections.defaultdict(lambda: 0)
    with open(args.bedtools_exons, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            rdID = line[RSIDCOL]
            in_exon[rdID] = 1
    
    #open bedtools overlap for whole transcript
    #store value for whether or not the SNP overlaps a reference gene    
    in_reference = collections.defaultdict(lambda: 0)
    with open(args.bedtools_trans_ref, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            rdID = line[RSIDCOL]
            in_reference[rdID] = 1
    
    #open bedtools overlap for whole transcript
    #store value for whether or not the SNP overlaps a reference gene    
    in_ref_exon = collections.defaultdict(lambda: 0)
    with open(args.bedtools_exons_ref, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            rdID = line[RSIDCOL]
            in_ref_exon[rdID] = 1
    
    #open bedtools overlap for compendia intergenic
    #store value for whether or not the SNP overlaps a compendia intergenic gene    
    in_intergenic = collections.defaultdict(lambda: 0)
    with open(args.bedtools_intergenic, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            rdID = line[RSIDCOL]
            in_intergenic[rdID] = 1
    
    
    #make a list of "cancer SNPs" 
    cancer_SNPs = set()
    with open(args.cancer_catalog, 'r') as f:
        headers = f.next().strip().split('\t')
        for row in f:
            row = row.strip().split('\t')
            if row != ['']:
                ID = row[headers.index('Snp_id_current')]
                if not ID: 
                    continue
                chrom = row[headers.index('Chr_id')]
                if not chrom: 
                    continue
                cancer_SNPs.add(ID)
    
    fileo = open("gwas_rsID_proc.tsv", 'w')
    with open(args.gwas_rsID_file, 'r') as f:
        headers = f.next().strip().split('\t')
        headero = headers + ['overlaps_compendia_gene', 
                             'overlaps_compendia_exon', 
                             'overlaps_ref_gene',
                             'overlaps_ref_exon',
                             'overlaps_compendia_intergenic',
                             "cancer_snp"]
        print >> fileo, '\t'.join(map(str, headero))
        for row in f:
            row = row.strip().split('\t')
            ID = row[headers.index('Snp_id_current')]
            in_comp = in_compendia[ID]
            in_ex = in_exon[ID]
            in_ref = in_reference[ID]
            in_ref_ex = in_ref_exon[ID]
            in_inter = in_intergenic[ID]
            if ID in cancer_SNPs: 
                cancer = 1
            else: 
                cancer = 0
            lineo = row + [in_comp, 
                           in_ex, 
                           in_ref,
                           in_ref_ex,
                           in_inter,
                           cancer]
            
            print >> fileo, '\t'.join(map(str, lineo))
    fileo.close()
            
    
            
    return 0



if __name__ == '__main__':

    sys.exit(main())