'''
Created on May 12, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np




def meta_read(line, header):
    metas = {}
    for item in header: 
        metas[item] = line[header.index(item)]
    return metas
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("meta")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.debug('reading metadata')
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_print = meta_header[:]
    meta_print.append('Lineage ES')
    meta_print.append('Lineage Percentile')
    meta_print.append('Lineage Type (Cancer/Normal)')
    meta_print.append('Cancer Vs Normal Es')
    meta_print.append('Cancer Vs Normal Percentile')
    meta_print.append('Average ES')
    meta_print.append('Average Percentile')
    
    print '\t'.join(meta_print)
    i = 0
    for line in meta_fh:
        i+=1
        if i%10000==0:
            logging.debug('Finished %d lines' % i) 
        line = line.strip().split('\t')
        metas = meta_read(line, meta_header)
        t_id = metas['transcript_id']
        func_type = metas['func_type']
        func_cat = metas['func_cat']
        func_dir = metas['func_dir']
        
        l_fam = 'NA'
        
        if func_cat == 'clat':
            if func_dir == 'up':
                l_es = metas['es.ctnt.up']
                l_frac = metas['frac.ctnt.up']
                cn_es = metas['es.cn.up']
                cn_frac = metas['frac.cn.up']
            if func_dir == 'dn':
                l_es = metas['es.ctnt.dn']
                l_frac = metas['frac.ctnt.dn']
                cn_es = metas['es.cn.dn']
                cn_frac = metas['frac.cn.dn']    
            a_es = str((float(l_es) + float(cn_es))/2.0)
            a_frac = str((float(l_frac) + float(cn_frac))/2.0)
        elif func_cat == 'cat':
            if func_dir == 'up':
                l_es = 'NA'
                l_frac = 'NA'
                cn_es = metas['es.cn.up']
                cn_frac = metas['frac.cn.up']
            if func_dir == 'dn':
                l_es = 'NA'
                l_frac = 'NA'
                cn_es = metas['es.cn.dn']
                cn_frac = metas['frac.cn.dn']    
            a_es = cn_es
            a_frac = cn_frac
        elif func_cat == 'at':
            if func_dir == 'up':
                l_es = metas['es.ctnt.up']
                l_frac = metas['frac.ctnt.up']
                l_fam = metas['ssea_set_family.ctnt.up']
                cn_es = 'NA'
                cn_frac = 'NA'
            if func_dir == 'dn':
                l_es = metas['es.ctnt.dn']
                l_frac = metas['frac.ctnt.dn']
                l_fam = metas['ssea_set_family.ctnt.dn']
                cn_es = 'NA'
                cn_frac = 'NA'    
            a_es = l_es
            a_frac = l_frac
        else: 
            l_es = 'NA'
            l_frac = 'NA'
            cn_es = 'NA'
            cn_frac = 'NA'
            a_es = 'NA'
            a_frac = 'NA'
        line.append(l_es)
        line.append(l_frac)
        line.append(l_fam)
        line.append(cn_es)
        line.append(cn_frac)
        line.append(a_es)
        line.append(a_frac)    
        print '\t'.join(line)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
