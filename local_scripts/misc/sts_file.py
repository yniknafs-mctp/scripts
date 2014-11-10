'''
Created on Apr 3, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np





    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("resp_ids")
    parser.add_argument("sts")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    resp_fh = open(args.resp_ids)
    resp_header = resp_fh.next().strip().split('\t')
    ids = []
    for line in resp_fh: 
        line = line.strip().split('\t')
        resp_id = line[resp_header.index('respondent_id')]
        ids.append(resp_id)
    
    sts_fh = open(args.sts)
    sts_headers = []
    for x in range(0,3):
        line = sts_fh.next().strip()
        sts_headers.append(line)
    temp = sts_fh.next().strip()
    
    #print header 
    for header in sts_headers: 
        print header + '\r'
    for id in ids: 
        print id.join(temp.split('####')) + '\r'
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
