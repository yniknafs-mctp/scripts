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
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_file")
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
        pt = line[header.index('tcga_legacy_id')]
        if pt.startswith("TCGA"):
            pt_split = pt.split('-')
            pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
            pt_to_RNA[pt_new] = RNA_sample
    
    
    #make dict to count the mutations for each gene 
    logging.debug('making KRAS mutant sample set')
    mut_dict = collections.defaultdict(lambda: [])
    fileh = open(args.maf_file)
    header = fileh.next().strip().split('\t')
    mut_pts = []
    for line in fileh: 
        line = line.strip().split('\t')
        func = line[header.index("Variant_Classification")]
        can_type = line[header.index("Cancer_Type")]
        aa_change = line[header.index("Protein_Change")]
        gene = line[header.index("Hugo_Symbol")]
        pt_split = line[header.index("Tumor_Sample_Barcode")].split('-')
        pt = '-'.join([pt_split[0], pt_split[1], pt_split[2], pt_split[3]])
        pt_rna = pt_to_RNA[pt]
        if pt_rna != "NOT_IN_COMPENDIA" and can_type == 'COADREAD':
            mut_pts.append(pt_rna)
        if func not in acceptable_muts: 
            continue
        if gene == 'KRAS' and can_type == 'COADREAD' and aa_change.startswith('p.G12'):
            mut_dict[gene].append(pt_rna)
    mut_pts = set(mut_pts)
    headero = ['set_name', 'set_description'] + list(mut_pts)
    print '\t'.join(headero)
    lineo = ['kras_luad_muts', 'LUAD pts with KRAS mutation']
    kras_muts = set(mut_dict['KRAS'])
    for x in xrange(2, len(headero)):
        if headero[x] in kras_muts:
            lineo.append(1)
        else: 
            lineo.append(0)
    print '\t'.join(map(str, lineo))
    logging.debug(len(mut_pts))
    logging.debug(len(mut_dict['KRAS']))
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
