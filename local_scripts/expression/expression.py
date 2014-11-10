'''
Created on Jan 20, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import glob

from ssea.lib.countdata import BigCountMatrix

'''
reads the memmap for the compendia expression data 
reports stats for the expression of each transcript 
stats: mean, max, min, 5th percentile, 95th percentile, median
'''

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("memmap_dir")
    parser.add_argument("metadata_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #make dict of lengths for count --> FPKM conversion    
    meta_file = open(args.metadata_file)
    meta_header = meta_file.next().strip().split('\t')
    t_meta_dict = {}
    logging.info("Parsing metadata file for transcript lengths")
    for line in meta_file: 
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        t_meta_dict[t_id] = line
    
        
    bcm = BigCountMatrix.open(args.memmap_dir)
    rows = bcm.rownames
    cols = bcm.colnames
    
    matrix = bcm.counts
    matrix_t = bcm.counts_t
    
    
    #loop through matrix by column to get total # of reads for FPKM conversion
    logging.info('Looping through matrix_t for number of reads per sample')
    tot_reads_array = []
    for x in xrange(len(matrix_t[:,1])):
        if (x%500) == 0: 
            logging.info('Finished %d/%d samples' % (x, len(matrix_t[:,1])))
        col_expr_array = matrix_t[x,]
        tot_reads_array.append(np.nansum(col_expr_array))
    tot_reads_np_array = np.array(tot_reads_array)
    
    #looping through matrix and converting to FPKM then reporting stats
    logging.info('Converting to FPKM then reporting stats')
    new_fields = [
                  'max',
                  '99.99th',
                  '99.9th',
                  '99.5th',
                  '99th',
                  '95th',
                  '90th',
                  'mean',
                  'median'
                  ]
    print '\t'.join(meta_header + new_fields)
    for x in xrange(len(rows)):
        if (x%5000) == 0: 
            logging.info('Finished %d/%d transcripts' % (x, len(rows))) 
        t_id = rows[x]
        metadata = t_meta_dict[t_id]
        t_len = metadata[meta_header.index('transcript_length')]
        t_len = int(t_len)
        count_array = matrix[x,]
        #convert to FPKM 
        fpkm_array = (count_array*10e8)/(t_len*tot_reads_np_array)
        max = fpkm_array.max()
        median = np.median(fpkm_array)
        mean = np.mean(fpkm_array)
        _99_99 = np.percentile(fpkm_array, 99.99)
        _99_9 = np.percentile(fpkm_array, 99.9)
        _99_5 = np.percentile(fpkm_array, 99.5)
        _99 = np.percentile(fpkm_array, 99)
        _95 = np.percentile(fpkm_array, 95)
        _90 = np.percentile(fpkm_array, 90)
        lineo = metadata + [max, _99_99, _99_9, _99_5, _99, _95, _90, mean, median]
        print '\t'.join(map(str,lineo))

#     
    return 0

if __name__ == '__main__': 
    sys.exit(main())
