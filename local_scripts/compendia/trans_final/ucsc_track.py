'''
Created on Mar 3, 2014
@author yniknafs
'''

'''
takes all named transcripts and produces a UCSC browser track 
with the functional names. The transcript is colored according 
to the class of transcript
'''


import os
import sys
import argparse
import logging
import numpy as np
import json
import collections

def meta_read(line, header):
    metas = {}
    for item in header: 
        metas[item] = line[header.index(item)]
    return metas

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("bed_file")
    parser.add_argument("meta_full")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    #print the bead file header
    print "track name=\"" + 'mitranscriptome' '\" description=\"' + 'mitranscriptome_named' + ' type=bedDetail'
        
    #read transcript metadata file to create dictionary for 
    logging.info('reading metadata')
    meta_fh = open(args.meta_full)
    meta_header = meta_fh.next().strip().split('\t')
    i=0
    t_set = []
    name_dict = {}
    cat_dict = {}
    type_dict = {}
    for line in meta_fh: 
        i+=1
        if i%10000==0:
            logging.debug('Finished %d lines' % i)
        line = line.strip().split('\t')
        metas = meta_read(line, meta_header)
        t_id = metas['transcript_id']
        t_set.append(t_id)
        func_name = metas['func_name']
        name_dict[t_id] = func_name
        func_cat = metas['func_cat']
        cat_dict[t_id] = func_cat
        func_type = metas['func_type']
        type_dict[t_id] = func_type
    t_set = set(t_set)
    
    
    result_file = open(args.bed_file)
    for line in result_file: 
        line = line.strip().split('\t')
        #inlcude line below if gene_id is in bed file
#         t_name = line[3].split('|')[1]
        t_id = line[3]
        if t_id not in t_set:
            continue
        name = name_dict[t_id]
        line[3] = name
        score = 1000
        line[4]=score
#         line.append(name)
#         line.append(desc)
        print '\t'.join(map(str,line))
            
    return 0

if __name__ == '__main__': 
    sys.exit(main())

