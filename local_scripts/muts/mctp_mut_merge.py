'''
Created on Mar 20, 2014
@author yniknafs
'''

'''
reads list of mctp patients and a key for the cancer type 
and concatenates all into a big mutation matrix
'''


import os
import sys
import argparse
import logging
import collections


headero = [
           'cancer_type',
           'patient_id',
           'hugo_symbol',
           'mutation_function',
           'codon_change',
           'protein_change'
           ]

EXON_IDX = 0
GENE_IDX = 1
FUNC_IDX = 2
AA_IDX = 3
CODON_IDX = 3
PROT_IDX = 4
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("mut_dir")
    parser.add_argument("mut_list")
    parser.add_argument("sample_key")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # print header
    print '\t'.join(headero)
    #read sample_info file to make dict for cohort each sample belongs to 
    fileh = open(args.sample_key)
    header = fileh.next().strip().replace('\"','').split(",")
    cohort_dict = collections.defaultdict(lambda: [])
    for line in fileh: 
        line = line.strip().replace('\"','').split(',')
        if len(line) < 20:
            continue
        cohort = line[header.index('Cancer Origin Site')]
        if cohort == '':
            continue
        pt_id = line[header.index("Patient or Cell Line ID")]
        cohort_dict[pt_id].append(cohort)
    
    cohort_dict_trim = collections.defaultdict(lambda: 'UNKNOWN')
    for key in cohort_dict.iterkeys():
        cohort_dict_trim[key] = ''.join(set(cohort_dict[key]))
    
    for mut_file in open(args.mut_list):
        mut_file = mut_file.strip()
        fileh = os.path.join(args.mut_dir, mut_file)
        pt = mut_file.split('.')[0]
        ca_type = cohort_dict_trim[pt]
        for line in open(fileh):
            line = line.strip().replace('\"','').split(',')
            exon = line[EXON_IDX]
            if exon != 'exonic':
                continue
            gene = line[GENE_IDX]
            func = line[FUNC_IDX]
            aa_long = line[AA_IDX].split(':')
            if len(aa_long) < 2: 
                continue
            codon = aa_long[CODON_IDX]
            prot = aa_long[PROT_IDX]
            lineo = [ca_type, 
                     pt,
                     gene,
                     func,
                     codon,
                     prot
                     ]
            
            print '\t'.join(map(str, lineo))
    
    
    
        
    
    
    
    
    
    logging.info('Running main script')
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
