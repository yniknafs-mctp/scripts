'''
Created on Mar 19, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections

'''
reads list of samples from pankaj
uses a csv file with sample info to 
group the samples
'''    

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("samples_list")
    parser.add_argument("samples_info")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    #read sample list file to get list of samples
    pts = []
    for sample in open(args.samples_list):
        pt = sample.split('.')[0]
        pts.append(pt)
    pts = set(pts)
    
    #read sample_info file to report the cohort each sample belongs to 
    fileh = open(args.samples_info)
    header = fileh.next().strip().replace('\"','').split(",")
    cohort_dict = collections.defaultdict(lambda: [])
    for line in fileh: 
        line = line.strip().replace('\"','').split(',')
        if len(line) < 20:
            continue
        cohort = line[header.index('Cancer Origin Site')]
        if cohort == '':
            continue
        pt_id = line[header.index("Patient or Cell Line ID")]
        cohort_dict[cohort].append(pt_id)
    
    for key in cohort_dict.iterkeys():
        count = len(set(cohort_dict[key]))
        print '\t'.join([key, str(count)])
        
        
    return 0

if __name__ == '__main__': 
    sys.exit(main())
