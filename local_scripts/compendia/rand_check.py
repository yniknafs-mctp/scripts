'''
Created on Apr 13, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import json



'''
investigate why bladder cancer versus normal is not working
'''



    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("json")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    for line in open(args.json):
        items = json.loads(line)
        seed = items['rand_seed']
        if seed == None:
            continue
        else: 
            seed = int(seed)
    ull
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
