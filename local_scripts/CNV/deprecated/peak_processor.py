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
    parser.add_argument("lesions_file")
    parser.add_argument("DNA_id_file")
    parser.add_argument("RNA_id_file")
    parser.add_argument('-a', dest="amp_threshold", default=0.1)
    parser.add_argument('-b', dest="del_threshold", default=0.1)

    args = parser.parse_args()
    
    a = os.path.abspath(args.lesions_file)
    b = b = os.path.dirname(a)
    c = os.path.basename(b)

    #building dictionary to convert GISTIC id to library ID   
    pt_to_RNA = collections.defaultdict(lambda: [])
    with open(args.RNA_id_file, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            RNA_sample, pt = row
            if pt.startswith("TCGA"):
                pt_split = pt.split('-')
                pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2]])
                line = [RNA_sample, pt_new]
                pt_to_RNA[pt_new] = RNA_sample
    DNA_to_RNA = collections.defaultdict(lambda: [])
    with open(args.DNA_id_file, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            DNA_sample, pt = row
            DNA_to_RNA[DNA_sample] = pt_to_RNA[pt][:]
            
    
    outfile_dir = b + "_lesions_processed.txt"
    print outfile_dir
    outfile = open(outfile_dir, 'w')
    
    with open(args.lesions_file, 'r') as f:
        header_fields = f.next().strip().split('\t') 
        reader=csv.reader(f,delimiter='\t')
        
        header_fields_copy = header_fields[:]
        header_fields_copy[9:9] = ["CNV_pos", "CNV_pos_c","CNV_null","CNV_null_c"]
        outfile.write('\t'.join(map(str, header_fields_copy)))
        outfile.write('\n')
        thresh_dict = collections.defaultdict(lambda: [])
        thresh_dict["Amplification"] = (args.amp_threshold, 1)
        thresh_dict["Deletion"] = (args.del_threshold, -1)
        for row in reader:
            CNV_null = []
            CNV_pos = []
            CNV_null_c = []
            CNV_pos_c = []
            row.pop()
            stripped_CNV = row[0].strip().replace(" ", "").split('-')
            stripped_ad = row[0].strip().split(' ')[0]
            if len(stripped_CNV) > 1:
                for x in xrange(9, len(row)):
                    product = float(row[x])*float(thresh_dict[stripped_ad][1])
                    if product > thresh_dict[stripped_ad][0]:
                        CNV_pos.append(header_fields[x])
                        if DNA_to_RNA[str(header_fields[x])] != []:
                            CNV_pos_c.append(DNA_to_RNA[str(header_fields[x])])
                    else: 
                        CNV_null.append(header_fields[x])
                        if DNA_to_RNA[str(header_fields[x])] != []:
                            CNV_null_c.append(DNA_to_RNA[str(header_fields[x])])
                row[9:9] = [CNV_pos, CNV_pos_c, CNV_null, CNV_null_c]
                outfile.write('\t'.join(map(str, row)))        
                outfile.write("\n")
        
    outfile.close()
#             
    
    
    logging.info("Finished")   

    return 0



if __name__ == '__main__':

    sys.exit(main())