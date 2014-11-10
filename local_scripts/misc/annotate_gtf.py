'''
Created on Apr 9, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
from assemblyline.lib.gtf import GTFFeature

'''
read in gtf and metadata
add fields in gtf for metadata
'''


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    parser.add_argument("meta_file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read metadata and make dict 
    logging.debug('reading metadata')
    meta_fh = open(args.meta_file)
    meta_header = meta_fh.next().strip().split("\t")
    meta_dict = {}
    for line in meta_fh: 
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        meta_dict[t_id] = line
    
    #loop through gtf and add meta fields
    logging.debug('appending metadata to gtf')
    for f in GTFFeature.parse(open(args.gtf_file)):
        t_id = f.attrs['transcript_id']
        meta = meta_dict[t_id]
        tcat = meta[meta_header.index('tcat')]
        tstatus = meta[meta_header.index('tstatus')]
        tgenic = meta[meta_header.index('tgenic')]
        f.attrs['tcat'] = tcat
        f.attrs['tstatus'] = tstatus
        f.attrs['tgenic'] = tgenic
        print str(f)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
