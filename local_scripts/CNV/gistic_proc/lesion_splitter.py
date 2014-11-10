'''
Created on October 22, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv

logging.basicConfig(level=logging.DEBUG, 
                    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

def read_lines(filename):
    lines = []
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    for line in fileh:
        if not line:
            continue
        if line.startswith("#"):
            continue
        line = line.strip()
        if not line:
            continue
        fields = line.strip().split('\t')
        lines.append(fields)
    return header_fields, lines


def process_lesion_file(lesion_file):
    
    logging.info("Processing lesion file: " + os.path.basename(lesion_file))
    lesions_header, lesions_raw = read_lines(lesion_file)
    
    #split lesionfile
    lesions = []
    
    for line in lesions_raw:
        peakID = line[lesions_header.index('Unique Name')]
        peakID_split = peakID.replace(' ','').split('-')
        peakID_clean = peakID_split[0]
        if len(peakID_split) > 1: 
            line[lesions_header.index('Unique Name')] = peakID_clean
            peak_data = peakID.split(' ')
            peak_type = peak_data[0]
            lesions.append(line)
            '''
                elif peak_type == "Deletion":
                    del_lesions.append(line)
                else: 
                    logging.info("Not valid peak type")
                    sys.exit()
            '''
    return lesions_header, lesions
