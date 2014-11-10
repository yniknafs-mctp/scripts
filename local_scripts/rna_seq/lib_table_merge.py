'''
Created on Jul 23, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np


'''
merges two library tables
'''


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("lib_table1")
    parser.add_argument("lib_table2")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Merging library tables')
    lib1_fh = open(args.lib_table1)
    lib1_dict = {}
    lib1_lib = []
    lib1_header = lib1_fh.next().strip().split('\t')
    for line in lib1_fh:
        line = line.strip().split('\t')
        lib_id = line[lib1_header.index('library_id')]
        lib1_lib.append(lib_id)
        lib1_dict[lib_id] = line
    lib1_lib = set(lib1_lib)
    
    lib2_fh = open(args.lib_table2)
    lib2_dict = {}
    lib2_lib = []
    lib2_header = lib2_fh.next().strip().split('\t')
    for line in lib2_fh:
        line = line.strip().split('\t')
        lib_id = line[lib2_header.index('library_id')]
        lib2_lib.append(lib_id)
        lib2_dict[lib_id] = line
    lib2_lib = set(lib2_lib)
    
    if lib1_header != lib2_header:
        logging.error("Both files given are not library tables")
        return 1
    
    dup_samples = lib1_lib & lib2_lib
    if len(dup_samples) == 0:
        print '\t'.join(lib1_header)
        for sample in lib1_dict.iterkeys():
            print '\t'.join(lib1_dict[sample])
        for sample in lib2_dict.iterkeys():
            print '\t'.join(lib2_dict[sample])
    
    dup_check = 0
    if len(dup_samples) != 0:
        logging.info('%d duplicates: %s' %(len(dup_samples), ','.join(list(dup_samples))))
        for sample in dup_samples:
            if lib1_dict[sample] != lib2_dict[sample]:
                dup_check+=1
                logging.error("Sample discrepancy for %s" % sample)
        if dup_check !=0: 
            return 1
        else: 
            print '\t'.join(lib1_header)
            for sample in lib1_dict.iterkeys():
                print '\t'.join(lib1_dict[sample])
            for sample in lib2_dict.iterkeys():
                if sample not in dup_samples:
                    print '\t'.join(lib2_dict[sample])
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
