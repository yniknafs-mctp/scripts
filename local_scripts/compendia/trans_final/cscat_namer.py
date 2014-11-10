'''
Created on Apr 15, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections


'''
reads in cscat meta file and appends a column with 
functional name (func_name)
'''

no_lat = [
          'bladder',
          'breast',
          'kich',
          'kirc',
          'kirp',
          'head_neck',
          'liver',
          'luad',
          'lusc',
          'stomach',
          'prostate',
          'thyroid',
          'kidney',
          'lung'
          ]

no_lat = set(no_lat)


def namer(type):
    if type == 'bladder':
        pre = 'BL'
    elif type == 'breast':
        pre = 'BR'
    elif type == 'kich':
        pre = 'KCH'
    elif type == 'kirc':
        pre = 'KC'
    elif type == 'kirp':
        pre = 'KP'
    elif type == 'head_neck':
        pre = 'HN'
    elif type == 'liver':
        pre = 'LV'
    elif type == 'luad':
        pre = 'LA'
    elif type == 'lusc':
        pre = 'LS'
    elif type == 'ovarian':
        pre = 'OV'
    elif type == 'stomach':
        pre = 'ST'
    elif type == 'prostate':
        pre = 'P'
    elif type == 'thyroid':
        pre = 'TH'
    elif type == 'aml': 
        pre = 'AM'
    elif type == 'cervical': 
        pre = 'CV'
    elif type == 'cml': 
        pre = 'CM'
    elif type == 'colorectal': 
        pre = 'CR'
    elif type == 'gbm': 
        pre = 'GB'
    elif type == 'heart': 
        pre = 'HR'
    elif type == 'hesc': 
        pre = 'ES'
    elif type == 'kidney': 
        pre = 'KD'
    elif type == 'lgg': 
        pre = 'LG'
    elif type == 'lung': 
        pre = 'LN'
    elif type == 'medulloblastoma': 
        pre = 'MB'
    elif type == 'melanoma': 
        pre = 'ME'
    elif type == 'mpn': 
        pre = 'MP'
    elif type == 'ovarian': 
        pre = 'OV'
    elif type == 'pancreatic': 
        pre = 'PN'
    elif type == 'skeletal_muscle': 
        pre = 'SM'
    elif type == 'uterine': 
        pre = 'UT'
    else: 
        print type
    return pre

def meta_read(line, header):
    metas = {}
    for item in header: 
        metas[item] = line[header.index(item)]
    return metas
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("meta")
    parser.add_argument("hugo")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.debug('reading hugo symbols')
    lncrna = []
    hugos = []
    hugo_fh = open(args.hugo)
    hugo_header = hugo_fh.next().strip().split('\t')
    for line in hugo_fh: 
        line = line.strip().split('\t')
        hugo = line[hugo_header.index('Approved Symbol')]      
        hugos.append(hugo)
        prev_hugo = line[hugo_header.index('Previous Symbols')].split(',')
        for symb in prev_hugo:
            hugos.append(symb)
        family = line[hugo_header.index('Gene Family Tag')].split(',')
        if 'LNCRNA' in family:
            lncrna.append(hugo.upper())
            

    logging.debug('reading metadata')
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_dict = {}
    suff_counter = collections.defaultdict(lambda: 1)
    clat_dict = collections.defaultdict(lambda: [])
    cat_dict = collections.defaultdict(lambda: [])
    lat_dict = collections.defaultdict(lambda: [])
    hiclincs = []
    clats = []
    name_dict = collections.defaultdict(lambda: 'NA')
    type_dict = collections.defaultdict(lambda: 'NA')
    category_dict = collections.defaultdict(lambda: 'NA')
    dir_dict = collections.defaultdict(lambda: 'NA')
    name_list = []
    nonsinglets = []
    i = 0
    ref_genes = []
    for line in meta_fh:
        i+=1
        if i%10000==0:
            logging.debug('Finished %d lines' % i) 
        line = line.strip().split('\t')
        metas = meta_read(line, meta_header)
        t_id = metas['transcript_id']
        meta_dict[t_id] = metas
        g_id = metas['gene_id']
        tcat = metas['tcat']
        #only use lncrna and tucp
        if tcat != 'lncrna' and tcat != 'tucp':
            continue
        ref_name = metas['ref_gene_name'].split(',')
        for name in ref_name: 
            ref_genes.append(name.upper())
        
        #parse data for all up-sets 
        ct_up_frac = metas['frac.ctnt.up']
        if ct_up_frac == 'NA': 
            ct_up_frac = 0
        ct_up_frac = float(ct_up_frac)
        ct_up_type = metas['ssea_type.ctnt.up']
        cn_up_frac = metas['frac.cn.up']
        if cn_up_frac == 'NA': 
            cn_up_frac = 0
        cn_up_frac = float(cn_up_frac)
        cn_up_type = metas['ssea_type.cn.up']       
        
        #parse data for all dn-sets 
        ct_dn_frac = metas['frac.ctnt.dn']
        if ct_dn_frac == 'NA':
            ct_dn_frac = 0
        ct_dn_frac = abs(float(ct_dn_frac))
        ct_dn_type = metas['ssea_type.ctnt.dn']
        cn_dn_frac = metas['frac.cn.dn']
        if cn_dn_frac == 'NA':
            cn_dn_frac = 0
        cn_dn_frac = abs(float(cn_dn_frac))
        cn_dn_type = metas['ssea_type.cn.dn']
        
        #identify clats
        clat_up = 0
        if ct_up_frac >= .95 and cn_up_frac >= .95:
            if ct_up_type == cn_up_type: 
                clat_up = (ct_up_frac + cn_up_frac)/2.0
                clat_up_type = ct_up_type
        clat_dn = 0
        if ct_dn_frac >= .95 and cn_dn_frac >= .95:
            if ct_dn_type == cn_dn_type: 
                clat_dn = (ct_dn_frac + cn_dn_frac)/2.0
                clat_dn_type = ct_dn_type
        
        clat_type = 'null'
        clat_val = 0
        if clat_up != 0 and clat_dn != 0:
            if clat_up >= clat_dn: 
                clat_type = clat_up_type
                clat_val = clat_up
            else: 
                clat_type = clat_dn_type
                clat_val = -1*clat_dn
        elif clat_up != 0 and clat_dn == 0: 
            clat_type = clat_up_type
            clat_val = clat_up
        elif clat_up == 0 and clat_dn != 0: 
            clat_type = clat_dn_type
            clat_val = -1*clat_dn
        
        #if a transcript does not make the clat cutoffs, check to see if it is a cat
        cat_type = 'null'
        cat_val = 0
        if clat_type == 'null':            
            if cn_up_frac >= .99 and cn_up_frac>cn_dn_frac: 
                cat_type = cn_up_type
                cat_val = cn_up_frac
            elif cn_dn_frac >= .99 and cn_dn_frac>cn_up_frac: 
                cat_type = cn_dn_type
                cat_val = -1*cn_dn_frac
        
        #if transcript does not make clat or cat cutoffs, check to see if it is a lat
        lat_type = 'null'
        lat_val = 0
        if cat_type == 'null' and clat_type == 'null':
            if ct_up_frac >= .99:
                lat_type = ct_up_type
                lat_val = ct_up_frac
            if ct_dn_frac >= .99:
                lat_type = ct_dn_type
                lat_val = -1*ct_dn_frac
        if lat_type in no_lat: 
            lat_type = 'null'
        
#         if t_id == 'T380286':
#             print 'what the cock'
#             print lat_type
#             print lat_val
#             print cat_val
#             print clat_val
        
        
        #finally if it does not make the clat or cat or lat cutoff, check to see if it has an UCE 
        hiclinc = 0
        uce =  metas['uce']
        if clat_type == 'null' and cat_type == 'null' and lat_type == 'null' and uce=='TRUE': 
            hiclinc = 1
        
        #store all named transcripts in respective dictionaries
        if clat_type != 'null': 
            value = [t_id, g_id, clat_val]
            cat_dict[clat_type].append(value)
            clats.append(t_id)
        if cat_type != 'null': 
            value = [t_id, g_id, cat_val]
            cat_dict[cat_type].append(value)
        if lat_type != 'null': 
            value = [t_id, g_id, lat_val]
            lat_dict[lat_type].append(value)
        if hiclinc == 1: 
            chrom = metas['chrom']
            if chrom == 'chrX': 
                chrom = 'chr26'
            if chrom == 'chrY': 
                chrom = 'chr26'
            start = metas['start']
            value = [t_id, g_id, int(chrom[3:]), int(start)]
            hiclincs.append(value)
    
    # make a set for all names already in use
    ref_genes = set(ref_genes)
    hugos = set(hugos)
    lncrna = set(lncrna)
    clats = set(clats)
    renamed = collections.defaultdict(lambda: 'null')
    #function to name the clats, cats, lats
    def type_namer(dict, label):
        for type in dict.iterkeys():
            if type == 'pancancer': 
                continue
            g_dict = collections.defaultdict(lambda: [])
            pre = namer(type)
            mid = label.upper()
            for rank, [t_id, g_id, score] in enumerate(sorted(dict[type], key = lambda x: (-x[2]))):
                g_dict[g_id].append([t_id, rank, score])
            g_rank = []
            for key in g_dict.iterkeys():
                rank = g_dict[key][0][1]
                g_rank.append([key, rank])
            for rank_final, [g_id, rank_rough] in enumerate(sorted(g_rank, key = lambda x: (x[1]))):
                trans = g_dict[g_id]
                annot_counter = 0
                for t_id,r, s in trans:
                    meta = meta_dict[t_id]
                    category = meta['category']
                    ref_name = set(meta['ref_gene_name'].upper().split(','))
                    if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
                        annot_counter +=1
                lnc_counter = len(trans) - annot_counter
                
                #make sure name that's generated is not already in the reference list
                suff_check = str(rank_final + 1)
                name_check = pre + mid + suff_check
                j=0
                while name_check in hugos: 
                    j +=1
                    suff_check = str(rank_final + 1 + j)
                    name_check = pre + mid + suff_check
                if lnc_counter > 0: 
                    hugos.add(name_check)
                
                for i,[t_id,r, s] in enumerate(trans): 
                    #check to see if it is an already named lncrna (from hugo) 
                    meta = meta_dict[t_id]
                    category = meta['category']
                    ref_name = set(meta['ref_gene_name'].upper().split(','))
                    if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
                        ref_name = list(ref_name.intersection(lncrna))[0]
                        name = ref_name + '.' + str(suff_counter[ref_name])
                        renamed[t_id] = '0'
                        if suff_counter[ref_name] == 2:
                            nonsinglets.append(ref_name)
                        suff_counter[ref_name] += 1
                    else: 
                        suff = suff_check + '.' + str(suff_counter[name_check])
                        if suff_counter[name_check] == 2:
                            nonsinglets.append(name_check)
                        suff_counter[name_check] +=1
                        name = pre + mid + suff
                        renamed[t_id] = '1'
                    if s >= 0:
                        dir_dict[t_id] = 'up'
                    else:
                        dir_dict[t_id] = 'dn'
                    name_dict[t_id] = name
                    type_dict[t_id] = type
                    if t_id in clats:
                        category_dict[t_id] = 'clat'
                    else: 
                        category_dict[t_id] = label.lower()
                    name_list.append(name)
                    
    
   
    #name the clats, cats, lats
    type_namer(clat_dict, 'cat')
    type_namer(cat_dict, 'cat')
    type_namer(lat_dict, 'at')
        
        
    #name the hiclincs
    g_dict = collections.defaultdict(lambda: [])
    mid = 'HICLINC'
    for rank, [t_id, g_id, chrom, start] in enumerate(sorted(hiclincs, key = lambda x: (x[2], x[3]))):
        g_dict[g_id].append([t_id, rank])
    g_rank = []
    for key in g_dict.iterkeys():
        rank = g_dict[key][0][1]
        g_rank.append([key, rank])
    for rank_final, [g_id, rank_rough] in enumerate(sorted(g_rank, key = lambda x: (x[1]))):
        trans = g_dict[g_id]
        annot_counter = 0
        for t_id,r in trans:
            meta = meta_dict[t_id]
            category = meta['category']
            ref_name = set(meta['ref_gene_name'].upper().split(','))
            if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
                annot_counter +=1
        lnc_counter = len(trans) - annot_counter
        
        #make sure name that's generated is not already in the reference list
        suff_check = str(rank_final + 1)
        name_check = mid + suff_check
        j=0
        while name_check in hugos: 
            j +=1
            suff_check = str(rank_final + 1 + j)
            name_check = mid + suff_check
        if lnc_counter > 0: 
            hugos.add(name_check)
        
        for i,[t_id,r] in enumerate(trans): 
            if t_id == 'T135389':
                continue
            #check to see if it is an already named lncrna (from hugo) 
            meta = meta_dict[t_id]
            category = meta['category']
            ref_name = set(meta['ref_gene_name'].upper().split(','))
            if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
                ref_name = list(ref_name.intersection(lncrna))[0]
                name = ref_name + '.' + str(suff_counter[ref_name])
                renamed[t_id] = '0'
                if suff_counter[ref_name] == 2:
                    nonsinglets.append(ref_name)
                suff_counter[ref_name] += 1
            else: 
                suff = suff_check + '.' + str(suff_counter[name_check])
                if suff_counter[name_check] == 2:
                    nonsinglets.append(name_check)
                suff_counter[name_check] +=1
                name = mid + suff
                renamed[t_id] = '1'
            name_dict[t_id] = name
            type_dict[t_id] = 'hiclinc'
            category_dict[t_id] = 'hiclinc'
            name_list.append(name)

    nonsinglets = set(nonsinglets)
    for key in name_dict.iterkeys():
        name = name_dict[key].split('.')[0]
        if name not in nonsinglets: 
            name_dict[key] = name
    
        
    logging.debug('printing new metadata')
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_header.append('func_name')
    meta_header.append('func_type')
    meta_header.append('func_cat')
    meta_header.append('func_dir')
    meta_header.append('renamed')
    print '\t'.join(meta_header)
    i = 0
    ref_genes = []
    for line in meta_fh:
        i+=1
        if i%10000==0:
            logging.debug('Printed %d lines' % i) 
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        func_name = name_dict[t_id]
        func_type = type_dict[t_id]
        func_cat = category_dict[t_id]
        func_dir = dir_dict[t_id]
        rename = renamed[t_id]
        line.append(func_name)
        line.append(func_type)
        line.append(func_cat)
        line.append(func_dir)
        line.append(rename )
        print '\t'.join(map(str, line))
    return 1
        
        
        
if __name__ == '__main__': 
    sys.exit(main())
