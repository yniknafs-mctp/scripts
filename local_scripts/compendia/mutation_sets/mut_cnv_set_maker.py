'''
Created on Mar 20, 2014
@author yniknafs
'''

'''
reads pancan maf file and pancan cnv file to generate sets for ssea
'''


import os
import sys
import argparse
import logging
import collections
import numpy as np

acceptable_muts = [                   
                    'Frame_Shift_Del',
                    'Frame_Shift_Ins',
                    'In_Frame_Del',
                    'In_Frame_Ins',
                    'Missense_Mutation',
                    'Nonsense_Mutation',
                    'Nonstop_Mutation',
                   ]

AMP_THRESH = 0.9
DEL_THRESH = -0.5
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_file")
    parser.add_argument("cnv_file")
    parser.add_argument("gene_list")
    parser.add_argument("pt_rna_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read lib file to convert TCGA pt ID to rna-seq lib id
    logging.info('Creating ID conversion dictionary')
    pt_to_RNA = collections.defaultdict(lambda: 'NOT_IN_COMPENDIA')
    fileh = open(args.pt_rna_file)
    header = fileh.next().strip().split("\t")
    for line in fileh: 
        line = line.strip().split('\t')
        RNA_sample = line[header.index('library_id')]
        pt = line[header.index('tcga_legacy_sample_id')]
        if pt.startswith("TCGA"):
            pt_split = pt.split('-')
            pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
            pt_to_RNA[pt_new] = RNA_sample
    
    #read gene list and make list
    logging.debug('reading gene list')
    gene_list = []
    for line in open(args.gene_list):
        line = line.strip()
        gene_list.append(line)
    gene_list = set(gene_list)
    
    #read cnv file and store arrays for the genes of interest 
    logging.debug('reading cnv file to store arrays for the genes of interest')
    fileh = open(args.cnv_file)
    cnv_pts_long = fileh.next().strip().split('\t')[3:]
    cnv_pts = []
    cnv_pts_set = []
    for pt in cnv_pts_long:
        pt_split = pt.split('-')
        pt = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
        pt_rna = pt_to_RNA[pt]
        cnv_pts.append(pt_rna)
        if pt_rna != []:
            cnv_pts_set.append(pt_rna)
    cnv_pts = np.array(cnv_pts)
    up_dict = collections.defaultdict(lambda: [])
    down_dict = collections.defaultdict(lambda: [])
    for line in fileh: 
        line = line.strip().split("\t")
        gene = line[0]
        if gene not in gene_list: 
            continue
        vals = np.array(line[3:]).astype(float)
        ups = vals > AMP_THRESH
        downs = vals < DEL_THRESH
        up_pts = cnv_pts[ups]
        down_pts = cnv_pts[downs]
        up_dict[gene] = set(up_pts)
        down_dict[gene] = set(down_pts)
    
    #make dict to count the mutations for each gene 
    logging.debug('counting mutations for each gene')
    mut_dict = collections.defaultdict(lambda: [])
    fileh = open(args.maf_file)
    header = fileh.next().strip().split('\t')
    mut_pts = []
    mut_pts_set = []
    for line in fileh: 
        line = line.strip().split('\t')
        func = line[header.index("Variant_Classification")]
        pt_split = line[header.index("Tumor_Sample_Barcode")].split('-')
        pt = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
        pt_rna = pt_to_RNA[pt]
        mut_pts.append(pt_rna)
        if pt_rna != []:
            mut_pts_set.append(pt_rna)
        if func not in acceptable_muts: 
            continue
        gene = line[header.index("Hugo_Symbol")]
        mut_dict[gene].append(pt_rna)
    mut_pts = set(mut_pts)
    
    
    #identify lists of pts for the combo sets
    logging.debug('Merging muts and cnvs')
    combo_up_dict = {}
    combo_down_dict = {}
    for gene in gene_list:
        pts = mut_dict[gene] #all pts with mut for this gene
        up_pts = up_dict[gene] #all pts with amp for this gene
        pt_muts_up = set(list(pts) + list(up_pts)) #combo of amp + mut
        down_pts = down_dict[gene] #all pts with del for this gene
        pt_muts_down = set(list(pts) + list(down_pts))#combo of del + mut
        combo_up_dict[gene] = pt_muts_up
        combo_down_dict[gene] = pt_muts_down
        

    with open('muts.smt', 'w') as f1:
        #print header for the muts file 
        mut_headero = ['set_name', 'set_description'] + list(mut_pts)
        print >>f1, '\t'.join(mut_headero)
        with open('muts_cnv.smt', 'w') as f2: 
            #print header for muts/cnv file
            mut_cnv_headero = ['set_name', 'set_description'] + list(set((mut_pts) & set(cnv_pts)))
            
            print >>f2, '\t'.join(mut_cnv_headero)
            i=0
            for gene in gene_list:
                i+=1
                logging.debug('Making sets for gene %d/%d: %s' %(i, len(gene_list),gene))
                #make set for just muts
                mut_name = 'pancan_' + gene + '_mutation'
                mut_desc = 'Pts with point mutation in %s' % gene
                lineo = [mut_name, mut_desc]
                for x in xrange(2, len(mut_headero)):
                    if mut_headero[x] in set(mut_dict[gene]):
                        lineo.append(1)
                    else:
                        lineo.append(0)
                print >>f1, '\t'.join(map(str, lineo))
                
                #make set for muts and amp
            
                mut_name = 'pancan_' + gene + '_mutation_amp'
                mut_desc = 'Pts with point mutation or amplification in %s' % gene
                lineo = [mut_name, mut_desc]
                for x in xrange(2,len(mut_cnv_headero)):
                    if mut_cnv_headero[x] in set(combo_up_dict[gene]):
                        lineo.append(1)
                    else:
                        lineo.append(0)
                print >>f2, '\t'.join(map(str, lineo))
                
                #make set for muts and del
                mut_name = 'pancan_' + gene + '_mutation_del'
                mut_desc = 'Pts with point mutation or deletion in %s' % gene
                lineo = [mut_name, mut_desc]
                for x in xrange(2, len(mut_cnv_headero)):
                    if mut_cnv_headero[x] in set(combo_down_dict[gene]):
                        lineo.append(1)
                    else:
                        lineo.append(0)
                print >>f2, '\t'.join(map(str, lineo))
                
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
