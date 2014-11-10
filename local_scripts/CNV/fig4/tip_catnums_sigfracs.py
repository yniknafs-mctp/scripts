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

'''
takes file containing all transcripts in peaks (w/ NES stats for that peak)
will collapse transcripts into genes, and output a file by peak
will output info about the number of transcripts in each peak that are 
of a particular category or are significant
'''
'''
input: tip_processed file from tip_ssea_metadata.py
THRESH = threshold used to identify "significant" genes
'''

THRESH = .001

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("transcripts_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #collapse into genes
    fileh = open(args.transcripts_file)
    header = fileh.next().strip().split('\t')
    gene_dict = collections.defaultdict(lambda: ['']*len(header))
    #saves every line if you ever want to do a per transcript analysis
    ###trans_dict = collections.defaultdict(lambda: [])
    
    for line in fileh:
        line = line.strip().split('\t')
        t_name = line[header.index('t_id')]
        peak_id = line[header.index('peak_id')]
        peak_type = line[header.index('peak_type')]
        gene = line[header.index('gene_id')]
        nes = line[header.index('nes')]
        nes = float(nes)
        #make key for dictionary that will store best gene line by peak
        key = (gene, peak_id)
        #store each line into the transcript dict
        ###key_t = (t_name, peak_id)
        ###trans_dict[key_t] = line
        if gene_dict[key][header.index('nes')] == '':
            gene_dict[key] = line
        if peak_type == "Amp":
            if nes > gene_dict[key][header.index('nes')]:
                gene_dict[key] = line
        if peak_type == "Del":
            if nes <= gene_dict[key][header.index('nes')]:
                gene_dict[key] = line
    #parse through gene_dict and make a dictionary with   
    #keys=peaks and values= all info for genes in that peak
    gip = collections.defaultdict(lambda: [])
    for g_name,peak in gene_dict.iterkeys():
        gip[peak].append(gene_dict[(g_name, peak)])
    
    #loop through gip dict (key = peak, value = genes in peak)
    #identify values for output file and print file
    fileo = open('peak_catnums_sigfracs.tsv', 'w')
    headers = ['peak_id', 
               'peak_type', 
               'best_nes', 
               'best_category', 
               'tot_genes',
               'sig_genes',
               'tot_prot_genes',
               'sig_prot_genes',
               'tot_lnc_genes',
               'sig_lnc_genes',
               'tot_nonprot_genes',
               'sig_nonprot_genes',
               'mean_nes_prot',
               'mean_nes_lnc',
               'mean_nes_nonprot',
               'mean_nes_prot_sig',
               'mean_nes_lnc_sig',
               'mean_nes_nonprot_sig',
               'mean_nes_prot_sig_sign',
               'mean_nes_lnc_sig_sign',
               'mean_nes_nonprot_sig_sign']  
    print >>fileo, '\t'.join(headers)
    for peak in gip.iterkeys():
        tot_genes = 0
        sig_genes = 0
        tot_prot_genes = 0
        sig_prot_genes = 0
        tot_lnc_genes = 0
        sig_lnc_genes = 0
        tot_nonprot_genes = 0
        sig_nonprot_genes = 0
        nes_list = []
        #list of all NES by category
        nes_list_nonprot = []
        nes_list_lnc = []
        nes_list_prot = []
        # list of NES for just sig genes
        nes_list_nonprot_sig = []
        nes_list_lnc_sig = []
        nes_list_prot_sig = []
        #following three enforce that NES is in direction of CNV
        nes_list_nonprot_sig_sign = []
        nes_list_lnc_sig_sign = []
        nes_list_prot_sig_sign = []
        cat_list = []
        for gene in gip[peak]:
            tot_genes+=1
            fdr = gene[header.index('fdr')]
            fdr = float(fdr)
            peak_type = gene[header.index('peak_type')] 
            nes = gene[header.index('nes')]
            nes = float(nes)
            cat = gene[header.index('category')]
            nes_list.append(nes)
            cat_list.append(cat)
            #count categories
            if cat == 'protein_coding':
                tot_prot_genes+=1
                nes_list_prot.append(nes)
            if cat == 'lncRNA':
                tot_lnc_genes+=1
                nes_list_lnc.append(nes)
            if cat != 'protein_coding': 
                tot_nonprot_genes+=1
                nes_list_nonprot.append(nes)
            #count sigs
            if fdr < THRESH: 
                if cat == 'protein_coding':
                    nes_list_prot_sig.append(nes)
                if cat == 'lncRNA':
                    nes_list_lnc_sig.append(nes)
                if cat != 'protein_coding':
                    nes_list_nonprot_sig.append(nes)
                if peak_type == "Amp":
                    if nes > 0:
                        sig_genes+=1
                        if cat == 'protein_coding':
                            sig_prot_genes+=1
                            nes_list_prot_sig_sign.append(nes)
                        if cat == 'lncRNA':
                            sig_lnc_genes+=1
                            nes_list_lnc_sig_sign.append(nes)
                        if cat != 'protein_coding': 
                            sig_nonprot_genes+=1
                            nes_list_nonprot_sig_sign.append(nes)
                if peak_type == "Del":
                    if nes < 0:
                        sig_genes+=1
                        if cat == 'protein_coding':
                            sig_prot_genes+=1
                            nes_list_prot_sig_sign.append(nes)
                        if cat == 'lncRNA':
                            sig_lnc_genes+=1
                            nes_list_lnc_sig_sign.append(nes)
                        if cat != 'protein_coding': 
                            sig_nonprot_genes+=1
                            nes_list_nonprot_sig_sign.append(nes)
        def means(list):
            if len(list) == 0: 
                mean = 0
            else: 
                mean = float(sum(list))/len(list)
            return mean
        mean_prot = means(nes_list_prot)
        mean_lnc = means(nes_list_lnc)
        mean_nonprot = means(nes_list_nonprot)
        mean_prot_sig = means(nes_list_prot_sig)
        mean_lnc_sig = means(nes_list_lnc_sig)
        mean_nonprot_sig = means(nes_list_nonprot_sig)
        mean_prot_sig_sign = means(nes_list_prot_sig_sign)
        mean_lnc_sig_sign = means(nes_list_lnc_sig_sign)
        mean_nonprot_sig_sign = means(nes_list_nonprot_sig_sign)
        
        
        
        #identify best gene category
        nes_list = np.array(nes_list)
        if peak_type == "Amp":
            best_index = nes_list.argmax()
        if peak_type == "Del":
            best_index = nes_list.argmin()
        
        
        best_cat = cat_list[best_index]
        best_nes = nes_list[best_index]
        
        lineo = [peak, peak_type, best_nes, best_cat, 
                 tot_genes, sig_genes, tot_prot_genes,
                 sig_prot_genes, tot_lnc_genes, 
                 sig_lnc_genes, tot_nonprot_genes, 
                 sig_nonprot_genes,
                 mean_prot, mean_lnc, mean_nonprot,
                 mean_prot_sig, mean_lnc_sig, mean_nonprot_sig,
                 mean_prot_sig_sign, mean_lnc_sig_sign, mean_nonprot_sig_sign]
        print >>fileo, '\t'.join(map(str,lineo))
    fileo.close()
    logging.info("Finished")
    return 0

if __name__ == '__main__': 
    sys.exit(main())
