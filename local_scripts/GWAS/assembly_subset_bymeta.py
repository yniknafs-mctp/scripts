'''
Created on Jan 16, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections

'''
takes the assembly metadata file and creates
a subset of the assembly bed file based off 
of categories of the metadata
'''



    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_bed")
    parser.add_argument("metadata_file")
    parser.add_argument("--fields", dest="fields",
                        help="comma separated list of fields to subset on (no spaces)")
    parser.add_argument("--values", dest="vals",
                        help='comma separated list of values to subset on (no spaces)' +
                        ' same order as fields')
    parser.add_argument("output_bed")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read in the fields and vals and make lists
    fields = args.fields.split(',')
    vals = args.vals.split(',')
    
    #check if fields and vals match up
    if len(fields) != len(vals):
        logging.info('Number of fields and vals given do not match...exiting')
        sys.exit()
    
    #check to make sure all the fields provided are in the metadata file
    meta = open(args.metadata_file)
    meta_header = meta.next().split('\t')
    meta_header_check = set(meta_header)
    fields_dict = collections.defaultdict(lambda: [])
    for x in xrange(len(fields)): 
        if fields[x] not in meta_header_check: 
            logging.info('Field \'%s\' not in metadata file, will not be used' % fields[x])
        else: 
            fields_dict[fields[x]].append(vals[x])
    #check to make sure at least one field given is in the metadata file
    if len(fields_dict.keys()) == 0: 
        logging.info('No provided fields match fields in metadata file...exiting')
        sys.exit()
    
    #print log for fields being used
    for field in fields_dict.iterkeys():
            val_options = set(fields_dict[field])
            logging.info('Using field \'%s\' with values \'%s\'' %(field, ','.join(val_options)))
    
    #store list of t_ids for the transcripts that will be subset from the bed file
    tIDs = set()
    for line in meta: 
        line = line.strip().split('\t')
        tID = line[meta_header.index('transcript_id')]
        check = set()
        for field in fields_dict.iterkeys():
            line_val = line[meta_header.index(field)]
            val_options = set(fields_dict[field])
            if line_val in val_options: 
                check.add(1)
            else: 
                check.add(0)
        if 0 not in check: 
            tIDs.add(tID)
    
    logging.info('%d transcripts match your criteria' % len(tIDs))
    
    #read BED file and print subset
    TIDCOL = 3
    bed = open(args.assembly_bed)
    fileo = open(args.output_bed, 'w')
    for line in bed: 
        line_split = line.strip().split('\t')
        tgID = line_split[TIDCOL]
        tID = tgID.split('|')[1]
        if tID in tIDs: 
            print >> fileo, line.strip()
    fileo.close()
        
    
    
    
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
