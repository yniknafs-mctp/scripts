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

from ssea.base import SampleSet, WeightVector
from ssea.config import Config
from ssea.algo import ssea_main
from ssea.report import Report


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

def timestamp():
    return datetime.fromtimestamp(time()).strftime('%Y-%m-%d-%H-%M-%S-%f')


def ssea_config(output_dir = "SSEA_%s" % (timestamp())):
    config = Config()
    # modify options as needed
    config.num_processors = 8
    config.name = 'myssea'
    config.weight_miss = 'log'
    config.weight_hit = 'log'
    config.weight_const = 1.1
    config.weight_noise = 0.05
    config.perms = 1000
    config.detailed_report_threshold = 0.05
    config.plot_conf_int = True
    config.conf_int = 0.95
    config.create_html = True
    config.create_plots = True
    config.output_dir = output_dir
    return config

def run_ssea(cnv_mat, null_mat, cnv_libs, null_libs, peakID, trans, output_dir):
    #prepare to run SSEA
        #must first combine the matrices for CNV and null    
    big_mat = np.hstack([cnv_mat, null_mat])
        #now combine the sampleID headers
    big_libs = cnv_libs + null_libs
    
    
    weight_vec_iter = WeightVector.from_data(rownames=trans,
                                     samples=big_libs,
                                     weight_matrix=big_mat, 
                                     na_value=-1)
    config = ssea_config(output_dir)
    
    # create a sample set
    sample_set = SampleSet(name=peakID,
                           desc='pts with CNV at ' + peakID,
                           value=set(cnv_libs))
    # pass list of sample sets to SSEA
    sample_sets = [sample_set]
    
    # run SSEA
    report_file = ssea_main(weight_vec_iter, sample_sets, config)
    
    out_list = []
    trans_list = []
    # parse output
    for rowdict in Report.parse(report_file):
        es = rowdict['es']
        nes = rowdict['nes']
        g_nes = rowdict['global_nes']
        global_qval = rowdict['global_fdr_q_value']
        available_fields = rowdict.keys()
        trans_name = rowdict['name']
        out_list.append([es, global_qval, nes, g_nes])
        trans_list.append(trans_name)
        
    return np.array(out_list).astype(float), trans_list


def ssea_ranker (ssea_list, q_val_cutoff, peak_type):
    if peak_type == "Amplification":
        sign = 1
    elif peak_type == "Deletion":
        sign = -1
    
    thresh=(ssea_list[:,1]<q_val_cutoff).astype(int)
    sign_adj = ssea_list*sign
    rank_big=np.lexsort((thresh, -sign_adj[:,0]))
    ranks_ranked = []
    for x in enumerate(rank_big):
        ranks_ranked.append(x)
    ranked = np.array(ranks_ranked)
    
    ranked_sorted = ranked[np.lexsort([ranked[:,0], ranked[:,1]])]
    ranks = []
    for x in ranked_sorted:
        ranks.append(x[0])
    
    
    ranks_out = []
    for x in xrange(len(ssea_list[:,0])):
        if thresh[x]==1: 
            ranks_out.append(ranks[x])
        else:
            ranks_out.append('not significant')
    return ranks_out

def ssea_dict(transcript_list, ssea_list, rank_list):
    dict = collections.defaultdict(lambda: [])
    for x in xrange(len(transcript_list)):
        dict[transcript_list[x]].append(rank_list[x])
        dict[transcript_list[x]].append(ssea_list[x][0])
        dict[transcript_list[x]].append(ssea_list[x][1])
        dict[transcript_list[x]].append(ssea_list[x][2])
        dict[transcript_list[x]].append(ssea_list[x][3])
    return dict
    
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
    parser.add_argument('-dt', dest="del_threshold", default=0.9)
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')
    
    #process the all_lesions file to separate header and remove redundant rows
    lesions_header, lesions = process_lesion_file(args.lesions_file)
    
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
    
    
    #open peak_file. This file will have information by peak 
    #(e.g. how many protein coding genes, ncrnas per peak)
    peak_file_dest = args.tumor_type + '_peaks.cnv'
    peak_file = open(peak_file_dest, 'w')
    
    #this file will have a peak name followed by two columns
    #one column will list the transcripts in the peak
    #the next will list the genes in the peak 
    items_peak_file_dest = args.tumor_type + '_items_in_peak.cnv'
    items_peak_file = open(items_peak_file_dest, 'w')
    
    #open transcript_file. This file will list information by transcript
    #It will include a rank for the transcript determined by ssea ES score
    transcript_file_dest = args.tumor_type + '_transcripts.cnv'
    transcript_file = open(transcript_file_dest, 'w')
    
    #add new fields for the peak_file
    new_lesion_header = lesions_header[:last_pheno_col]
    new_lesion_header.extend(['#of proteins',
                              '#of antisense',
                              '#of ncrna',
                              '#of intronic',
                              '#of intergenic'])
    
    #create header for the items_in_peak file
    items_header = ['peakID',
                    'transcripts in peak',
                    'genes in peak',
                    'CNV pts',
                    'null pts']

    #add new fields for the transcript_file                    
    transcript_header = ['cancer type', 
                         'peakID',
                         'peak_type',
                         'cnv_count',
                         'null_count',
                         'es',
                         'global FDR',
                         'nes',
                         'global_nes',
                         'rank_in_peak']     
    
    
    #write headers for peak_file and items_file to file
    print >>peak_file, '\t'.join(map(str,new_lesion_header)) 
    print >>items_peak_file, '\t'.join(map(str,items_header))
    
    i = 0
    
    #will loop through each line of the all_lesions file
    #each line is information for each peak
    for peak in lesions:
        #print trans_dict
        peakID = peak[0]    
        tran_inpeak = trans_dict[peakID]
        gene_inpeak = genes_dict[peakID]
        
        #returns a list of the sampleIDs for pts that do or do not have the CNV
        cnv_pts, null_pts = CNV_pts(lesions_header, peak, last_pheno_col, args.amp_threshold, args.del_threshold)
        #returns expression data for the CNV pts
        #cnv_mat has expression data
        #cnv_meta contains the meta_data for each transcript
        #trans returns a list of the transcripts with order preserved to the matrix
        #cnv_libs returns a list of the libIDs with order preserved to the matrix
        #meta_header returns the header_fields for the meta_data (used to print file)
        cnv_mat, cnv_meta, trans, cnv_libs, meta_header = get_expression(cnv_pts, tran_inpeak)
        
        #returns expression data for the null pts
        null_mat, _, _, null_libs, _ = get_expression(null_pts, tran_inpeak)
        
        #prints header line for the transcript_file. Need the meta_header to do so
        #so this must be performed within this loop, but only once
        if i==0:
            new_header = meta_header+transcript_header 
            print>>transcript_file, '\t'.join(map(str, new_header))
        i+=1
        
        #dictionary to count the transcript types in each peak
        #this information will get printed to the peak_file
        type_count_dict = {'protein':0, 'antisense':0, 'ncrna':0, 'intronic':0, 'intergenic':0}
        
        #count the transcript types for this peak
        for line in cnv_meta:
            trans_type = line[transcript_type_col]
            type_count_dict[trans_type]+=1
        
        #prepare line to be printed to the peak_file
        line_out = peak[:last_pheno_col]
        line_out.extend([type_count_dict["protein"],
                         type_count_dict["antisense"],
                         type_count_dict["ncrna"],
                         type_count_dict["intronic"],
                         type_count_dict["intergenic"]])
        
        #prepare line to be printed to the items_file
        line_out_items = [peakID, 
                          ','.join(tran_inpeak),
                          ','.join(gene_inpeak),
                          cnv_pts,
                          null_pts]
        
        #print line to the peak_file and items_file
        print >>peak_file, '\t'.join(map(str, line_out))
        print >>items_peak_file, '\t'.join(map(str,line_out_items))
        
        #define variables that will be used in the transcript_file
        #that do not change from transcript to transcript within the peak
        cancer_type = args.tumor_type
        peak_type = peakID.split('Peak')[0]
        cnv_count = len(cnv_libs)
        null_count = len(null_libs)
        
        #run SSEA
        ssea_list, trans_list = run_ssea(cnv_mat, null_mat, cnv_libs, null_libs, peakID, trans, peakID)
        #ssea ranker
        rank_list = ssea_ranker(ssea_list, args.qval, peak_type)
        #make ssea dict
        ssea_dictionary = ssea_dict(trans_list, ssea_list, rank_list)
         
        
        
        
        #print ssea_list
        #print is_sig
        #loop through each transcript in the peak
        
        for x in xrange(len(trans)):
            transcript = trans[x] 
            if ssea_dictionary[transcript] != []:
                #define variables to print in transcript_file that do change
                #from transcript to transcript
                
                rank = ssea_dictionary[transcript][rank_in]
                es = ssea_dictionary[transcript][es_in]
                qval = ssea_dictionary[transcript][fdr_in]
                nes = ssea_dictionary[transcript][nes_in]
                g_nes = ssea_dictionary[transcript][g_nes_in]
                trans_meta = cnv_meta[x]
                cnv_expr = cnv_mat[x, :]
                null_expr = null_mat[x,:]
                
                #prepare data to be printed to trnscript_flie
                transcript_info = [cancer_type, 
                                   peakID, 
                                   peak_type,
                                   cnv_count ,
                                   null_count, 
                                   es,
                                   qval,
                                   nes,
                                   g_nes,
                                   rank]
                
                line_trans = trans_meta + transcript_info
                
                #print line of transcript_file
                print >>transcript_file, '\t'.join(map(str,line_trans))
            
            
    peak_file.close()
    items_peak_file.close()
    transcript_file.close()
    

    
    
    logging.info('Finished')
    return 0

if __name__ == '__main__': 
    sys.exit(main())
