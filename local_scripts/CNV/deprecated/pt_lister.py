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
    logging.info("Running code")
    parser = argparse.ArgumentParser()
    parser.add_argument("processed_file")
    args = parser.parse_args()
    
    with open(args.processed_file, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        a = next(f).strip().split("\t") 
        newlist = a[10][:]
        newlist2 = newlist.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
        for x in newlist2: 
            print x
        newlist3 = a[12][:]
        newlist4 = newlist3.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
        for x in newlist4:
            print x
#             
    
    
    logging.info("Finished")   

    return 0



if __name__ == '__main__':

    sys.exit(main())