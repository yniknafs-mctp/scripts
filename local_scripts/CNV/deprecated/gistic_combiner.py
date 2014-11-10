import argparse
import logging
import os
import sys
import csv
import time
import collections
import bisect
import operator
import subprocess
import glob
import numpy
from scipy import stats
import plot_test
    
def main():

    
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    
    parser.add_argument("file_list")
    args = parser.parse_args()
    
    
    
    outfile = open('pan_can_seg.txt', 'w')
        
    for line in open(args.file_list):
        if not line:
            continue
        if line.startswith("#"):
            continue
        line = line.strip()
        if not line:
            continue
        for row in open(line):
            if not row:
                continue
            if row.startswith("#"):
                continue
            row = row.strip()
            if not row:
                continue
            outfile.write(row)
            outfile.write('\n')
    outfile.close()
        
    
        
       


    
    
    logging.info("Finished")   
    


    return 0



if __name__ == '__main__':

    sys.exit(main())