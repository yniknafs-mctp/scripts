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
    parser.add_argument("RNA_sampleID_file")
    parser.add_argument("DNA_sampleID_file")
    args = parser.parse_args()
    
    outfile_dir = args.RNA_sampleID_file[:-4] + "_trimmed.txt"
    print outfile_dir
    outfile = open(outfile_dir, 'w')
    
    pt_to_RNA = collections.defaultdict(lambda: [])
    #creating dictionary
    with open(args.RNA_sampleID_file, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            RNA_sample, pt = row
            if pt.startswith("TCGA"):
                pt_split = pt.split('-')
                pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2]])
                line = [RNA_sample, pt_new]
                pt_to_RNA[pt_new] = RNA_sample
                outfile.write('\t'.join(map(str, line)))
                outfile.write("\n")
    outfile.close()
    DNA_to_RNA = collections.defaultdict(lambda: [])
    with open(args.DNA_sampleID_file, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            DNA_sample, pt = row
            DNA_to_RNA[DNA_sample] = pt_to_RNA[pt][:]
        
        print DNA_to_RNA["08e895e6-2f45-436c-9b33-9768b66dc33b"]
            
    
    
    logging.info("Finished")   

    return 0



if __name__ == '__main__':

    sys.exit(main())