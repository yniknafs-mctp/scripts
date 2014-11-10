'''
Created on Mar 28, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections


'''
will take a list of genes and a list of patients
also takes a big MAF matrix and a CNV matrix 
also takes an expression matrix for a given gene/genes 
will output a file where each row is a patient and the columns are 
all the stats for muts, cnv, and expression
pt_list should be a list of TCGA-IDs
'''

CNV_LAST_PHENO = 3

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("pt_list")
    parser.add_argument("gene_list")
    parser.add_argument("maf_file")
    parser.add_argument("cnv_file")
    parser.add_argument("schlap")
    parser.add_argument("spink")
    parser.add_argument("id_conv")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #make a list of genes to query
    gene_list = []
    for line in open(args.gene_list): 
        line = line.strip()
        gene_list.append(line)
    
    # make dict for ID conversion
    #read lib file to convert TCGA pt ID to rna-seq lib id
    logging.info('Creating ID conversion dictionary')
    RNA_to_pt = collections.defaultdict(lambda: 'NOT_IN_COMPENDIA')
    fileh = open(args.id_conv)
    header = fileh.next().strip().split("\t")
    for line in fileh: 
        line = line.strip().split('\t')
        RNA_sample = line[header.index('library_id')]
        pt = line[header.index('tcga_legacy_id')]
        if pt.startswith("TCGA"):
            pt_split = pt.split('-')
            pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
            RNA_to_pt[RNA_sample] = pt_new
    
    
    logging.info('Reading expression data')
    schlap = {}
    spink = {}
    #read in schlap and spink expresion, store in dict
    schlap_fh = open(args.schlap)
    schlap_header = schlap_fh.next().strip().split('\t')
    schlap_expr = schlap_fh.next().strip().split('\t')
    for x in xrange(1, len(schlap_header)):
        pt_rna = schlap_header[x]
        pt_tcga = RNA_to_pt[pt_rna]
        pt_tcga = '-'.join(pt_tcga.split('-')[0:4])
        expr = schlap_expr[x]
        schlap[pt_tcga] = expr
    spink_fh = open(args.spink)
    spink_header = spink_fh.next().strip().split('\t')
    spink_expr = spink_fh.next().strip().split('\t')
    for x in xrange(1, len(spink_header)):
        pt_rna = spink_header[x]
        pt_tcga = RNA_to_pt[pt_rna]
        pt_tcga = '-'.join(pt_tcga.split('-')[0:4])
        expr = spink_expr[x]
        spink[pt_tcga] = expr
    
    # make a set for all patients being analyzed
    pt_set = []
    for pt in open(args.pt_list):
        pt = pt.strip()
        pt = '-'.join(pt.split('-')[0:4])
        pt_set.append(pt)
    pt_set = set(pt_set)
    
    # make a set for all genes being analyzed
    gene_list = []
    for gene in open(args.gene_list):
        gene = gene.strip()
        gene_list.append(gene)
    gene_set = set(gene_list)
    
    logging.info('Parsing mutation file')
    #read mut file to make dict of muts 
    mut_fh = open(args.maf_file)
    mut_header = mut_fh.next().strip().split("\t")
    mut_dict = collections.defaultdict(lambda: [])
    for line in mut_fh: 
        line = line.strip().split('\t')
        aa = line[mut_header.index("Protein_Change")]
        gene = line[mut_header.index("Hugo_Symbol")]
        pt = line[mut_header.index("Tumor_Sample_Barcode")]
        pt = '-'.join(pt.split('-')[0:4])
        if pt not in pt_set:
            continue
        if gene not in gene_set:
            continue
        key = (pt, gene)
        mut_dict[key].append(aa)
    
    logging.info('Parsing copy number file')
    #read cnv file to make a dict of cnv scores
    cnv_fh = open(args.cnv_file)
    cnv_header = cnv_fh.next().strip().split('\t')
    cnv_dict = collections.defaultdict(lambda: 'no_cnv_data')
    for line in cnv_fh:
        line = line.strip().split('\t')
        gene = line[cnv_header.index("Gene Symbol")]
        if gene not in gene_set:
            continue
        for x in xrange(CNV_LAST_PHENO, len(cnv_header)):
            pt = cnv_header[x]
            pt = '-'.join(pt.split('-')[0:4])
            if pt not in pt_set:
                continue 
            cn = float(line[x])
            if cn > 0.9: 
                value = 'amp('+ str(cn) + ')'
            elif cn < -0.3:
                value = 'del('+ str(cn) + ')'
            else: 
                value = ''
            key = (pt, gene)
            cnv_dict[key] = value
    
    #write output 
    headero = [
               'patient',
               'schlap_expression_fpkm',
               'spink_status'
               ]
    for gene in gene_list:
        headero.append(gene+'_cnv')
        headero.append(gene+'_mut')
    
    print '\t'.join(headero)
    
    for pt in pt_set:
        pt = '-'.join(pt.split('-')[0:4])
        lineo = [pt,schlap[pt],spink[pt]] 
        for gene in gene_list: 
            cnv = cnv_dict[(pt,gene)]
            mut = ','.join(mut_dict[(pt,gene)])
            lineo.append(cnv)
            lineo.append(mut)
        print '\t'.join(map(str, lineo))
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
