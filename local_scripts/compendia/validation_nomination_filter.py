'''
Created on Aug 5, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import collections



def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("validation_list")
    parser.add_argument("glist")
    parser.add_argument("sequences")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Reading g_list')
    glist = set()
    fileh = open(args.glist)
    header = fileh.next().strip().split('\t')
    for gid in fileh:
        gid = gid.strip()
        glist.add(gid)
        
    logging.info('Reading sequence data')
    seq_dict = collections.defaultdict(lambda: 'na')
    fileh = open(args.sequences)
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0].split('|')[1]
        seq = line[2].upper().replace('|','[]')
        seq_dict[tid] = seq
    
    fileh = open(args.validation_list)
    header = fileh.next().strip().split('\t')
    logging.debug(header)
    print '\t'.join(header)
    for line in fileh:
        line = line.strip().split('\t')
        fields = []
        for item in line: 
            if item != '':
                fields.append(item)
        fields = fields[:-1]
        if len(fields) <= 1: 
            continue
        tid = fields[header.index("transcript_id")]
        gid = fields[header.index("gene_id")]
        val = fields[9]
        if val == '0':
            if gid not in glist: 
                continue
        seq = seq_dict[tid]
        fields.append(seq)
        print '\t'.join(fields)
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
