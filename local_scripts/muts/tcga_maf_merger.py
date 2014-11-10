'''
Created on Mar 20, 2014
@author yniknafs
'''

'''
takes many MAFs from multiple centers and combines into one big maf
reads in a file that is a list of the MAF files in one column, and the types in another
'''


import os
import sys
import argparse
import logging
import numpy as np


header_common = [
               'Hugo_Symbol',
               'Center',
               'NCBI_Build',
               'Chromosome',
               'Start_position',
               'End_position',
               'Strand',
               'Variant_Classification',
               'Variant_Type',
               'Tumor_Sample_Barcode',
               'Verification_Status',
               'Verification_Status',
               'Mutation_Status'
               ]

header_out = [
               'Cancer_Type',
               'Hugo_Symbol',
               'Center',
               'NCBI_Build',
               'Chromosome',
               'Start_position',
               'End_position',
               'Strand',
               'Variant_Classification',
               'Variant_Type',
               'Codon_Change',
               'Protein_Change',
               'Tumor_Sample_Barcode',
               'Verification_Status',
               'Verification_Status',
               'Mutation_Status'
               ]

print '\t'.join(header_out)

def broad_change(header, line):
    codon = line[header.index('cDNA_Change')]
    prot = line[header.index('Protein_Change')]
    return codon, prot

def wust_change(header, line):
    for x in xrange(len(header)): 
        if header[x].endswith('_WU'): 
            header[x] = header[x][:-3]
    ref = line[header.index('reference')]
    var = line[header.index('variant')]
    c_pos = line[header.index('c_position')]
    codon = c_pos + ref + '>' + var
    prot = line[header.index('amino_acid_change')]
    return codon, prot

def bcm_change(header, line):
    prot = line[header.index('AAChange')]
    codon = line[header.index('ChromChange')]
    return codon, prot
            
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    #make sure maf_list has a second column for the ca_type
    parser.add_argument("maf_list")
    parser.add_argument("maf_dir")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    for line in open(args.maf_list):
        line = line.strip().split('\t')
        if len(line) < 2: 
            logging.error('Input file must specify cancer type for each maf file')
            sys.exit()
        fileh, ca_type = line 
        fileh = os.path.join(args.maf_dir, fileh)
        filer = open(fileh)
        for line in filer: 
            if line.startswith('#'):
                continue
            if line.startswith("Hugo"): 
                header = line.strip().split('\t')
                header_lower = line.strip().lower().split('\t')
                continue
            line = line.strip().split("\t")
            center = line[header.index('Center')]
            if center == 'broad.mit.edu':
                codon, prot = broad_change(header, line)
            elif center == 'genome.wustl.edu':
                codon, prot = wust_change(header, line)
            elif center == 'hgsc.bcm.edu':
                codon, prot = bcm_change(header, line)
            lineo = []
            for field in header_common: 
                lineo.append(line[header_lower.index(field.lower())])
            codon_index = header.index("Variant_Type")
            prot_index = codon_index + 1
            lineo.insert(codon_index, codon)
            lineo.insert(prot_index, prot)
            lineo = [ca_type] + lineo
            print '\t'.join(map(str, lineo))
    
    
    logging.info('Running main script')
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
