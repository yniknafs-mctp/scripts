'''
Created on Aug 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging





    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("lib_table")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Converting library table')
    headero = [
              'sample_id',
              'library_id',
              'gtf_file',
              'bam_file'
              ]
    print '\t'.join(headero)
    fileh = open(args.lib_table)
    header = fileh.next().strip().split('\t')
    for line in fileh:
        line = line.strip().split('\t')
        sample_id = line[header.index('sample_id')]
        library_id = line[header.index('library_id')]
        gtf_file = ''
        results_dir = line[header.index('results_dir')]
        bam_file = os.path.join(results_dir, 'tophat', 'accepted_hits.bam')
        if not os.path.exists(bam_file):
            logging.error("BAM File not Found: %d" % bam_file)
            continue
        lineo = [sample_id, library_id, gtf_file, bam_file]
        print '\t'.join(lineo)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

