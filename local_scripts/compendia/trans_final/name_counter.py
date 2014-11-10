'''
Created on Sep 11, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import collections




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("assoc_file")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    fileh = open(args.assoc_file)
    header = fileh.next().strip().split('\t')
    i = set()
    k = set()
    gene_dict = collections.defaultdict(lambda: set())
    l_dict = {}
    for line in fileh: 
        line = line.strip().split('\t')
        gid = line[header.index('gene_id')]
        ct_up = line[header.index('ct_up')]
        cn_up = line[header.index('cn_up')]
        clat_up = line[header.index('clat_up')]
        ct_dn = line[header.index('ct_dn')]
        cn_dn = line[header.index('cn_dn')]
        clat_dn = line[header.index('clat_dn')]
        k.add(gid)
        j = set()
        l = [ct_up, cn_up, clat_up, ct_dn, cn_dn, clat_dn]
        z = set()
        for item in l:                
            if item == 'NA':
                continue
            else: 
                for tissue in item.split(','):
                    gene_dict[gid].add(tissue)
                    z.add(tissue)
    
    for gene in gene_dict.iterkeys():
        j = gene_dict[gene]
        if len(j)==1:
            i.add(gene)
    print len(i)
    print len(k)
        
        
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
