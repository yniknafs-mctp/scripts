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

def main():
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("files_file")
    args = parser.parse_args()
    logging.info("Starting")
    
    with open(args.files_file, 'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            lesions_file, matrix_file, tumor_type = row
            p = subprocess.call(['python','/home/yniknafs/workspace/CNV/gistic_processor.py', lesions_file, matrix_file, '.', tumor_type])
            logging.info("Finished " + tumor_type)

    return 0



if __name__ == '__main__':

    sys.exit(main())