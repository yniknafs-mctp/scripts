'''
Created on Sep 26, 2013

@author: mkiyer
'''
'''
Created on May 1, 2013

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import numpy as np

EXPRESSION_DIR = "/mctp/users/mkiyer/projects/assemblyline/compendia_2013-May-13/gene_expression_matrix"

def read_lines(filename):
    lines = []
    for line in open(filename):
        if not line:
            continue
        if line.startswith("#"):
            continue
        line = line.strip()
        if not line:
            continue
        lines.append(line)
    return lines

def read_table(filename, primary_key='tracking_id', subset=None):
    if subset is None:
        subset = set()
    metadata = []
    inds = []
    keys = []
    fileh = open(filename)
    header_fields = fileh.next().strip().split('\t')
    id_ind = header_fields.index(primary_key)
    ind = 0
    for line in fileh:
        fields = line.strip().split('\t')
        key = fields[id_ind]        
        if len(subset) > 0 and key in subset:
            metadata.append(fields)
            keys.append(key)
            inds.append(ind)
        ind += 1
    fileh.close()
    return ind, header_fields, metadata, inds, keys

def get_expression(library_ids, trans_ids):
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    if not os.path.exists(EXPRESSION_DIR):
        parser.error("Expression matrix directory '%s' not found" % (args.matrix_dir))
    #logging.info("Getting peak expression")
    # read input library ids
    library_id_set = set(library_ids)
    # read input gene ids
    gene_id_set = set(trans_ids)

    # load gene expression
    pheno_file = os.path.join(EXPRESSION_DIR, 'phenos.txt')
    metadata_file = os.path.join(EXPRESSION_DIR, 'metadata.txt')
    matrix_file = os.path.join(EXPRESSION_DIR, 'isoform_fpkm.mmap')
    # find libraries
    ncols, lib_header_fields, lib_metadata, lib_inds, lib_ids = \
        read_table(pheno_file, primary_key='library_id', 
                   subset=library_id_set)
    if len(lib_inds) == 0:
        print "No libraries found"
        return 1
    # find genes
    nrows, g_header_fields, g_metadata, g_inds, g_ids = \
        read_table(metadata_file, primary_key='tracking_id', 
                   subset=gene_id_set)
    if len(g_inds) == 0:
        print "No genes found"
        return 1
    # read gene expression
    mat = np.memmap(matrix_file, dtype='float32', mode='r', 
                    shape=(nrows,ncols))
    # get subset of matrix
    submat = mat[g_inds,:]
    submat = submat[:,lib_inds]
    submat = submat.astype('float')
    submat[submat<0] = -1
    
    return submat, g_metadata, g_ids, lib_ids, g_header_fields
    

if __name__ == '__main__': 
    sys.exit(main())
