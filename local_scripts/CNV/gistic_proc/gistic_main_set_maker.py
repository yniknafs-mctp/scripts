'''
Created on October 28, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv
from lesion_splitter import process_lesion_file
from bed_tools import bed_handler
import subprocess
from cnv_to_seq_ID_conv import id_converter
from CNV_samples import CNV_pts
from get_isoform_expression_CNV import get_expression
from cnv_plotter import plotter
import collections
from datetime import datetime
from time import time
# import classes needed to set up SSEA




#index for last column of metadata for the all_lesions file
last_pheno_col = 9 

#index for the column in the compendia metadata file that has the transcript type information
transcript_type_col = 8

#index for rank in ssea_dict, es score, and global FDR qvalue
rank_in = 0
es_in = 1
fdr_in = 2
nes_in = 3
g_nes_in = 4


def lesion_namer(headers, lesions, type):
    #dictionary that will have the chr loc as a key and the values will be a list of peaks with that chr loc
    locs_dict = collections.defaultdict(lambda: [])
    
    #dictionary that will have the old peakID as a key and the new unique ID as a value
    new_name_dict = collections.defaultdict(lambda: [])
    
    #loop through the lesions file to pair up peaks with the chr loc
    for row in lesions: 
        peakID = row[headers.index('Unique Name')]
        peak_type = peakID.split('Peak')[0]
        locus = row[headers.index('Descriptor')].strip()
        tuple = (peak_type, locus)
        locs_dict[tuple].append(peakID)

    
    #loop through the chr locs to assign unique names to the peaks
    for key in locs_dict.iterkeys():
        peak_type, loc = key

        #list of peaks with the same chr loc
        peaks = locs_dict[key]
        type_abb = peak_type[0:3]
        
        
        if len(peaks) == 1:
            new_name = [type, type_abb, loc]
            new_name_dict[peaks[0]] = '_'.join(new_name)
            
        elif len(peaks) > 1: 
            for x, peak in enumerate(peaks):
                x+=1
                new_name = [type, type_abb, loc, str(x)]
                new_name_dict[peak] = '_'.join(new_name)
                
    return new_name_dict   
        
    
        
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("lesions_file")
    parser.add_argument("matrix_file")
    parser.add_argument("tumor_type")
    parser.add_argument("output_dir")
    parser.add_argument('-qval', dest="qval", default=0.05)
    parser.add_argument("-bed", dest="bed_file", default = '/home/yniknafs/misc_files/compendia_2013_hg18.bed')
    parser.add_argument('-at', dest="amp_threshold", default=1.2)
    parser.add_argument('-dt', dest="del_threshold", default=0.5)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    
    
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)


    
    
    #process the all_lesions file to separate header and remove redundant rows
    lesions_header, lesions = process_lesion_file(args.lesions_file)
    unique_dict = lesion_namer(lesions_header, lesions, args.tumor_type)
    
    #run bedtools intersect on all peaks to return a dictionary 
    #with each peak as a key and the transcripts/genes in that peak
    #as the value
    trans_dict, genes_dict = bed_handler(lesions_header, lesions, args.bed_file)
    
    #produce a dictionary to convert SNP6 sample IDs into rna_seq sampleIDs
    #that are found in the compendia
    id_conv_dict = id_converter(args.matrix_file)
    
    #convert the header fields of the all_lesions file from SNP6 to rnaseq IDs
    #note: this will leave '[]' for samples that do not have an rnaseq ID
    for x in xrange(last_pheno_col, len(lesions_header)):
        lesions_header[x] = id_conv_dict[lesions_header[x]]
    
    
    #open peak_file. This file will sets for all the peaks 
    peak_file_dest = os.path.join(args.output_dir,  args.tumor_type + '_peaks.smt')
    peak_file = open(peak_file_dest, 'w')
    
    #open transcript_file. This file will have transcripts and the peak they are in 
    transcript_file_dest = os.path.join(args.output_dir,  args.tumor_type + '_transcripts.tsv')
    transcript_file = open(transcript_file_dest, 'w')
    
    #open file for patient ids  for the universe 
    universe_dest = os.path.join(args.output_dir, args.tumor_type + '_ptIDs.txt')
    universe_file = open(universe_dest, 'w')
    
    cnv_pts, null_pts = CNV_pts(lesions_header, lesions[0], last_pheno_col, float(args.amp_threshold), float(args.del_threshold))
    universe = cnv_pts + null_pts
    
    #write pts in unverse to universe file
    for pt in universe: 
        print >> universe_file, pt
    universe_file.close()
    
    
    #will loop through each line of the all_lesions file
    #each line is information for each peak
    
    out_line = ['Set name', 'Set description'] + cnv_pts + null_pts
    peak_file.write('\t'.join(map(str,out_line))+'\n')
    for peak in lesions:
        peakID = peak[0]    
        upid = unique_dict[peakID]
        upid = upid.replace('_', '-')
        type = upid.split('-')[1]
        coord = peak[lesions_header.index('Wide Peak Limits')].split('(')[0]
        desc = 'Pts with ' + type + ' at ' + coord
        #returns a list of the sampleIDs for pts that do or do not have the CNV
        cnv_pts, null_pts = CNV_pts(lesions_header, peak, last_pheno_col, args.amp_threshold, args.del_threshold)
        cnv_pts = set(cnv_pts)
        membership = []
        for pt in universe: 
            if pt in cnv_pts: 
                membership.append(1)
            else: 
                membership.append(0)
        out_line = [upid, desc] + membership
        print >> peak_file, '\t'.join(map(str, out_line)) 
        

        #all transcripts in peak 
        tran_inpeak = trans_dict[peakID]
        for transcript in tran_inpeak:
            line = [transcript, upid]
            print >> transcript_file, '\t'.join(map(str, line))
    peak_file.close()
    transcript_file.close()

    logging.info('Finished')
    return 0

if __name__ == '__main__': 
    sys.exit(main())
