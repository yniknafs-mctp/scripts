'''
Created on Jul 23, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import subprocess
import collections

'''
Take a transcript metadata file, and add fields for expression stats for 
various cohorts
'''

from ssea.lib.countdata import BigCountMatrix

def parse_library_table(filename, names):
    namedict = dict((name,i) for i,name in enumerate(names))
    col_inds = []
    with open(filename) as f:
        header_fields = f.next().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            d = dict(zip(header_fields, fields))
            if d['library_id'] not in namedict:
                logging.warning('Library %s not found in matrix' % (d['library_id']))
                continue
            i = namedict[d['library_id']]
            col_inds.append(i)
    return col_inds

def parse_transcript_table(filename, names):
    namedict = dict((name,i) for i,name in enumerate(names))
    row_inds = []
    lengths = []
    with open(filename) as f:
        header_fields = f.next().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            d = dict(zip(header_fields, fields))
            t_id = d['transcript_id']
            if t_id not in namedict:
                logging.warning('Transcript %s not in matrix' % (t_id))
                continue
            i = namedict[t_id]
            length = int(d['transcript_length']) / 1000.0
            row_inds.append(i)
            lengths.append(length)
    return row_inds, lengths

    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_dir")
    parser.add_argument("title")
    parser.add_argument("library_ids")
    parser.add_argument("metadata")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    # open matrix
    bm = BigCountMatrix.open(args.matrix_dir)
    # get pheno data
    logging.debug('Reading library metadata')
    col_inds = parse_library_table(args.library_ids, bm.colnames)
    # get transcript metadata
    logging.debug('Reading transcript metadata')
    row_inds, lengths = parse_transcript_table(args.metadata, bm.rownames)
    
    bm.count_tots /= 1.0e6
    
    mean_dict = collections.defaultdict(lambda: float('NaN'))
    nf_dict = collections.defaultdict(lambda: float('NaN'))
    nn_dict = collections.defaultdict(lambda: float('NaN'))
    mean_dict_f = collections.defaultdict(lambda: float('NaN'))
    nf_dict_f = collections.defaultdict(lambda: float('NaN'))
    nn_dict_f = collections.defaultdict(lambda: float('NaN'))
    
    logging.debug('Calculating expression stats')
    k = 0 
    for i,lengthkb in zip(row_inds,lengths):
        k+=1
        if k%10000 == 0: 
            logging.debug('Finished %d transcripts' % k)
        
        t_id = bm.rownames[i]
        counts = bm.counts[i, col_inds]
        # ignore nans
        valid_inds = np.isfinite(counts)
        # normalize counts
        count_vals = (counts[valid_inds] / bm.size_factors[valid_inds])
        fpkm_vals = (counts[valid_inds] / (lengthkb*bm.count_tots[valid_inds]))
        
        mean = np.mean(count_vals)
        mean_f = np.mean(fpkm_vals)
        if list(count_vals) == []: 
            continue
        nf = np.percentile(count_vals, 95)
        nf_f = np.percentile(fpkm_vals, 95)
        nn = np.percentile(count_vals, 99)
        nn_f = np.percentile(fpkm_vals, 99)
        
        mean_dict[t_id] = mean
        mean_dict_f[t_id] = mean_f
        nf_dict[t_id] = nf
        nf_dict_f[t_id] = nf_f
        nn_dict[t_id] = nn
        nn_dict_f[t_id] = nn_f
    
    logging.debug('Printing output')
    meta_fh = open(args.metadata)
    meta_header = meta_fh.next().strip().split('\t')
    meta_header.append(args.title+'_count_mean')
    meta_header.append(args.title+'_count_95')
    meta_header.append(args.title+'_count_99')
    meta_header.append(args.title+'_fpkm_mean')
    meta_header.append(args.title+'_fpkm_95')
    meta_header.append(args.title+'_fpkm_99')
    print '\t'.join(meta_header)
    for line in meta_fh:
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        mean = mean_dict[t_id]
        mean_f = mean_dict_f[t_id]
        nf = nf_dict[t_id]
        nf_f = nf_dict_f[t_id]
        nn = nn_dict[t_id]
        nn_f = nn_dict_f[t_id]
        if mean == 'na' or mean_f == 'na':
            line.append('na')
        line.append("%.3f" % mean)
        line.append("%.3f" % nf)
        line.append("%.3f" % nn)
        line.append("%.3f" % mean_f)
        line.append("%.3f" % nf_f)
        line.append("%.3f" % nn_f)
        print '\t'.join(line)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
