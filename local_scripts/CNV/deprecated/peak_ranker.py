'''
Created on October 22, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv
import collections

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("cnv_file")
    parser.add_argument('peak_type')
    parser.add_argument("-sorter", dest='key', default='nes')
    parser.add_argument("-p", dest="p_cutoff", default=0.05)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    if args.peak_type == 'amp':
        sign = 1
    elif args.peak_type == 'del':
        sign = -1
    else:
        logging.info('Not a valid peak type')
        sys.exit()
    
    
    key_dict = collections.defaultdict(lambda: [])
    rank_dict = collections.defaultdict(lambda: [])
    with open(args.cnv_file, 'r') as f:
        header_fields = f.next().strip().split('\t')
        reader = csv.reader(f, delimiter='\t')
        for row in reader: 
            transcript = row[header_fields.index('transcript_id')]
            peak = row[header_fields.index('peakID')]
            key = row[header_fields.index(args.key)]
            pval = row[header_fields.index(args.key+'_pval')]
            if float(pval)<args.p_cutoff:
                key_dict[peak].append([key, transcript])
            else:
                rank_dict[transcript] = 'not_significant'
    
    for key in key_dict.iterkeys():
        trans_in_peak = key_dict[key]
        tip_array = np.array(trans_in_peak)
        sort_index = (tip_array).argsort(axis=0)
        sorted = tip_array[sort_index[:,0][::sign], :]
        for rank, item in list(enumerate(sorted)):
            rank +=1
            rank_dict[item[1]] = rank
        
    outfile_dest = os.path.join(os.path.dirname(args.cnv_file), 
                                'ranked_' + os.path.basename(args.cnv_file))
    outfile = open(outfile_dest, 'w')
    with open(args.cnv_file, 'r') as f:
        header_fields = f.next().strip().split('\t')
        add_loc = header_fields.index('log_nes_pvalplot_hyperlink')

        header_new = header_fields[:]
        header_new[add_loc:add_loc] = ['rank_in_peak']
        outfile.write(header_new + '\n')
        
        reader = csv.reader(f, delimiter='\t')
        for row in reader: 
            transcript=row[header_fields.index('transcript_id')]
            row[add_loc:add_loc]=rank_dict[transcript]
            outfile
            
    
    
    
    
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
