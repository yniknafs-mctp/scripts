'''
Created on Jan 18, 2014
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
    parser.add_argument("lesions_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Running main script')
    print args.lesions_file
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
