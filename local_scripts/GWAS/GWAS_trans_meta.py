'''
Created on Jan 20, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections

#BED columns
TRANS_CHR = 0
TRANS_START = 1
TRANS_STOP = 2
TRANS_ID = 3
SNP_CHR = 12
SNP_START = 13
SNP_STOP = 14
SNP_NAME = 15
SNP_DIST = 16


'''
takes BEDTools closest file that lists nearest snp to all of the 
compendia transcripts

also takes the exon intersect file

finally takes assembly metadata file 

appends columns for nearest snp, dist to nearest snp, and if the overlap is exonic
'''

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("closest_file")
    parser.add_argument("metadata_file")
    parser.add_argument("exon_intersect_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    snp_dict = collections.defaultdict(lambda: [])
    for line in open(args.closest_file):
        line = line.strip().split('\t')
        trans_chr = line[TRANS_CHR]
        trans_start = line[TRANS_START]
        trans_stop = line[TRANS_STOP]
        trans_id = line[TRANS_ID]
        snp_chr = line[SNP_CHR]
        if snp_chr == 'chr23':
            print 'tit'
            snp_chr = 'chrX'
        snp_start = line[SNP_START]
        snp_stop = line[SNP_STOP]
        snp_name = line[SNP_NAME]
        if snp_name != '.':
            snp_name = 'rs'+ snp_name
        snp_dist = line[SNP_DIST] 
        snp_dict[trans_id].append((snp_name, snp_dist)) 
    
    exon_dict = collections.defaultdict(lambda: 0)
    for line in open(args.exon_intersect_file):
        line = line.strip().split('\t')
        t_id = line[TRANS_ID]
        exon_dict[t_id] = 1
    
    meta_file = open(args.metadata_file)
    header = meta_file.next().strip().split('\t')
    header = header + ['snp_id', 'snp_dist', 'exonic_overlap']
    print '\t'.join(header) 
    for line in meta_file: 
        line = line.strip().split('\t')
        t_id = line[header.index('transcript_id')]
        snp_names = []
        snp_dists = []
        for snp_name, snp_dist in snp_dict[t_id]: 
            snp_names.append(snp_name)
            snp_dists.append(snp_dist)
        exon = exon_dict[t_id]        
        line = line + [','.join(snp_names), snp_dists[0], str(exon)]
        print '\t'.join(line)
        

    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
