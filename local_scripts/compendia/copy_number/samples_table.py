'''
Created on Mar 14, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np


'''
Goes through all GISTIC analyses and builds table of all samples
used, and adds a column if we had matched RNAseq data
'''


    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("files_file")
    parser.add_argument("rna_id_key")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    

    #build dictionary to key the rna library IDs to the patient ID
    key_f = open(args.rna_id_key)
    key_h = key_f.next().strip().split('\t')
    rna_id_dict = {}
    ca_type_dict = {}
    for line in key_f:
        line = line.strip().split('\t')
        if len(line) < 2:
            continue
        sample_id = line[key_h.index("library_id")]
        pt_id_long = line[key_h.index("description")]
#         ca_type = line[key_h.index("tcga_cancer_type")]
        pt_id_split = pt_id_long.split('-')
        if len(pt_id_split) < 3:
            continue
        pt_id = '-'.join([pt_id_split[0], pt_id_split[1], pt_id_split[2]])
        rna_id_dict[sample_id.lower()] = pt_id
#         ca_type_dict[pt_id] = ca_type
    
    #print header 
    header = [
                  'patient_id',
                  'sample_id',
                  'cancer_type',
                  'matched_rna_seq'
                  ]
    print '\t'.join(header)

    #read files file and parse through each cancer type 
    for line in open(args.files_file):
        line = line.strip().split('\t')
        lesions = line[0]
        matrix = line[1]
        catype = line[2]
        smt = line[3] 
        logging.info('Parsing ca_type %s' % catype)
        
        #build dictionary to key the cnv library IDs to the patient ID
        matrix_f = open(matrix)
        matrix_h = matrix_f.next().strip().split('\t')
        cnv_id_dict = {}
        for line in matrix_f:
            line = line.strip().split('\t')
            sample_id = line[matrix_h.index("sampleID")]
            pt_id = line[matrix_h.index("_PATIENT")]
            cnv_id_dict[sample_id.lower()] = pt_id
        
        if not os.path.exists(smt):
            logging.info('No RNA-seq for ca_type %s, SSEA not run' % catype)
        else: 
            #parse through smt file to make list of pt_ids for which we have matched RNA-seq
            seq_set = []
            smt_f = open(smt)
            smt_h = smt_f.next().strip().split('\t')
            #index for where metadata stops
            index = smt_h.index('Set description') + 1
            for lib_id in smt_h[index:]:
    #             if lib_id.lower() not in id_dict.keys(): 
    #                 print catype
    #                 continue
                pt_id = rna_id_dict[lib_id.lower()]
                seq_set.append(pt_id)
            seq_set = set(seq_set)
        
        #parse through all_lesions file to gather all samples used for copy number analysis
        lesions_f = open(lesions)
        lesions_h = lesions_f.next().strip().split('\t')
        #index for where metadata stops
        index = lesions_h.index('Amplitude Threshold') + 1
        for lib_id in lesions_h[index:]:
            pt_id = cnv_id_dict[lib_id.lower()]
#             catype = ca_type_dict[pt_id]
            
            if pt_id in seq_set: 
                seq = 1
            else: 
                seq = 0
            lineo = [pt_id, lib_id, catype, str(seq)]
            print '\t'.join(lineo)
            
            
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
