'''
Created on Nov 10, 2014
@author yniknafs
'''

'''
take metadata file and expression data. filter for 
candidate hiclincs for CRISPR screen
'''


import os
import sys
import argparse
import logging
import collections





    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

