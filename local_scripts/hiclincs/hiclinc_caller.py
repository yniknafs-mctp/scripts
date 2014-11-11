'''
Created on Nov 10, 2014
@author yniknafs
'''

'''
take metadata file that includes expression data. 
filter for candidate hiclincs for CRISPR screen
'''

'''
ssea file = /mctp/projects/mitranscriptome/rebuttal/finished_files/all.trans.ssea.txt
'''

import os
import sys
import argparse
import logging
import collections



FPKM_THRESH = 1.5



    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata")
    parser.add_argument("--ssea", dest = 'ssea', 
                        default = '/mctp/projects/mitranscriptome/rebuttal/finished_files/all.trans.ssea.txt')
    parser.add_argument("--cell-line", dest = 'cell_line',
                        default = '/mctp/users/yniknafs/projects/hiclincs/hiclinc.cell_line.expr.txt')
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    hiclincs = []
    
    cell_line_dict = collections.defaultdict(lambda: 'na')
    fileh = open(args.cell_line)
    cl_header = fileh.next().strip().split('\t')
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0]
        cell_line_dict[tid] = line[1:]
        hiclincs.append(tid)
        
    hiclincs = set(hiclincs)
    
    cn_up_dict = collections.defaultdict(lambda: 'na')
    cn_dn_dict = collections.defaultdict(lambda: 'na')
    ct_up_dict = collections.defaultdict(lambda: 'na')
    ct_dn_dict = collections.defaultdict(lambda: 'na')
    fileh = open(args.ssea)
    header = fileh.next().strip().split('\t')
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[0]
        cn_up_dict[tid] = line[header.index('cn_up')]
        cn_dn_dict[tid] = line[header.index('cn_dn')]
        ct_up_dict[tid] = line[header.index('ct_up')]
        ct_dn_dict[tid] = line[header.index('ct_dn')]
    
    
    fileh = open(args.metadata)
    header = fileh.next().strip().split('\t')
    headero = ['tid', 'gid', 'chrom', 'start', 'end', 'strand', 
               'cn_up', 'cn_dn', 'ct_up', 'ct_dn'] + cl_header[1:]
    print '\t'.join(headero)
    genes = []
    cn = []
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[0]
        gid = line[header.index('gene_id')]
        chrom = line[1]
        start = line[2]
        end = line[3]
        strand = line[4]
        expr_counter = 0 
        if tid not in hiclincs: continue
        for head in header: 
            if head.endswith('fpkm_99'): 
                expr = line[header.index(head)]
                if float(expr) > FPKM_THRESH: 
                    expr_counter+=1
        if expr_counter == 0: continue
        cn_up = cn_up_dict[tid]
        cn_dn = cn_dn_dict[tid]
        ct_up = ct_up_dict[tid]
        ct_dn = ct_dn_dict[tid]
        
        if cn_up.lower() != 'na' or cn_dn.lower() != 'na' or \
            ct_up.lower() != 'na' or ct_dn.lower() != 'na': cn.append(gid) 
            
        
         
        
        cl = cell_line_dict[tid]
        
        lineo = [tid, gid, chrom, start, end, strand,
                 cn_up, cn_dn, ct_up, ct_dn] + cl
        print '\t'.join(lineo)
        
        genes.append(gid)
    logging.debug(len(set(genes)))
    logging.debug(len(set(cn)))
        
    return 1
         




        
            
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

