'''
Created on Nov 18, 2013

@author: yniknafs
'''

import sys
import argparse
import logging
import collections
import gzip


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("CNV_file")
    parser.add_argument("thresh")
    parser.add_argument("reflat")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    file = gzip.open(args.reflat)
    chr_dict = collections.defaultdict(lambda: '')
    for line in file: 
        line = line.strip().split('\t')
        gene = line[0]
        chr = line[2]
        pq = int(line[4])
        pq_call = ''
        if pq < 45600000:
            pq_call = 'p'
        else: 
            pq_call = 'q'
        chr_dict[gene] = (chr, pq_call)
    
    
    file = open(args.CNV_file)
    header = file.next().strip().split('\t')
    tot = len(header) - 1
    thresh = float(args.thresh)
    for line in file: 
        line = line.strip().split('\t')
        gene = line[header.index('Gene Symbol')]
        del_count = 0
        for sample in line[1:]:
            cnv = float(sample)
            if cnv < thresh: 
                del_count +=1
        frac = del_count/float(tot)
        chr_pq = chr_dict[gene]
        if chr_pq:
            chr = chr_pq[0]
            pq = chr_pq[1]
        else: 
            chr = 'no meta'
            pq = 'no meta'
        lineo = [gene, chr, pq, del_count, frac]
        print '\t'.join(map(str,lineo)) 
    
    
        

if __name__ == '__main__':
    sys.exit(main())