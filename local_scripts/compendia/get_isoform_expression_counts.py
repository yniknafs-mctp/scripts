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

EXPRESSION_DIR = "/mctp/projects/ssea/isoform_count_matrix_v7"
from ssea.lib.countdata import BigCountMatrix


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

def read_table(filename, subset=None):
    if subset is None:
        subset = set()
    inds = []
    keys = []
    fileh = open(filename)
    ind = 0
    for line in fileh:
        key = line.strip()
        if len(subset) > 0 and key in subset:
            keys.append(key)
            inds.append(ind)
        ind += 1
    fileh.close()
    return ind, inds, keys

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("library_ids")
    parser.add_argument("transcript_ids")
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    if not os.path.exists(EXPRESSION_DIR):
        parser.error("Expression matrix directory '%s' not found" % (args.matrix_dir))
    # read input library ids
    library_ids = read_lines(args.library_ids)
    library_id_set = set(library_ids)
    # read input gene ids
    gene_ids = read_lines(args.transcript_ids)
    gene_id_set = set(gene_ids)

    # load gene expression   
    pheno_file = os.path.join(EXPRESSION_DIR, 'colnames.txt')
    metadata_file = os.path.join(EXPRESSION_DIR, 'rownames.txt')
    
    # find libraries
    ncols, lib_inds, lib_ids = \
        read_table(pheno_file, 
                   subset=library_id_set)
    if len(lib_inds) == 0:
        print "No libraries found"
        return 1
    # find genes
    logging.info('Acquiring data for %s libraries' % len(lib_ids))
    nrows, g_inds, g_ids = \
        read_table(metadata_file, 
                   subset=gene_id_set)
    if len(g_inds) == 0:
        print "No genes found"
        return 1
    logging.info('Acquiring data for %s transcripts' % len(g_ids))
    
    # read gene expression
    bm = BigCountMatrix.open(EXPRESSION_DIR)
    mat = bm.counts
    
    # get subset of matrix
    logging.info("performing gene subset")
    submat = mat[g_inds,:]
    logging.info("performing library subset")
    submat = submat[:,lib_inds]
    
    # write expr file
    output_expr_file = args.output_prefix + ".expr.tsv"
    fileh = open(output_expr_file, 'w')
    fields = ['gene_id']
    fields.extend(lib_ids)
    logging.debug('Printing expression file')
    print >>fileh, '\t'.join(fields)
    for i in xrange(len(g_inds)):
        if i%1000==0: 
            logging.debug("Finished %d/%d transcripts" % (i, len(g_inds)))
        fields = [g_ids[i]]
        for x in submat[i, :]:
            if x<0:
                fields.append('NA')
            else: 
                fields.append(str(x))
        print >>fileh, '\t'.join(fields)
    fileh.close()
    return 0

if __name__ == '__main__': 
    sys.exit(main())
