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
from operator import itemgetter

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
        pre = 'BLC'
    elif type == 'breast':
        pre = 'BRC'
    elif type == 'kich':
        pre = 'KHC'
    elif type == 'kirc':
        pre = 'KCC'
    elif type == 'kirp':
        pre = 'KPC'
    elif type == 'head_neck':
        pre = 'HNC'
    elif type == 'liver':
        pre = 'LVC'
    elif type == 'luad':
        pre = 'LAC'
    elif type == 'lusc':
        pre = 'LSC'
    elif type == 'ovarian':
        pre = 'OV'
    elif type == 'stomach':
        pre = 'STC'
    elif type == 'prostate':
        pre = 'PC'
    elif type == 'thyroid':
        pre = 'THC'
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
    elif type == 'embryonic_stem_cells': 
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
    parser.add_argument("assoc")
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
    hugos = set(hugos)
    
    logging.debug('reading metadata')
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_dict = {}
    hiclincs = []
    i = 0
    for line in meta_fh:
        i+=1
        if i%10000==0:
            logging.debug('Finished %d lines' % i) 
        line = line.strip().split('\t')
        metas = meta_read(line, meta_header)
        t_id = metas['transcript_id']
        meta_dict[t_id] = metas
        uce =  metas['uce']
        tcat = metas['tcat']
        if uce == 'TRUE' and tcat in ['lncrna', 'tucp']:
            hiclincs.append(t_id)
        
    
    logging.debug('reading associations file')
    assoc_gene_dict = collections.defaultdict(lambda: set())
    fileh = open(args.assoc)
    cats = set()
    header = fileh.next().strip().split('\t')
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[header.index('transcript_id')]
        gid = line[header.index('gene_id')]
        ct_frac_up = line[header.index('ct_frac_up')]
        cn_frac_up = line[header.index('cn_frac_up')]
        clat_frac_up = line[header.index('clat_frac_up')]
        ct_frac_dn = line[header.index('ct_frac_dn')]
        cn_frac_dn = line[header.index('cn_frac_dn')]
        clat_frac_dn = line[header.index('clat_frac_dn')]
        l = [ct_frac_up, cn_frac_up, clat_frac_up, ct_frac_dn, cn_frac_dn, clat_frac_dn]
#         if ct_frac_up == 'NA' and cn_frac_up == 'NA' and clat_frac_up == 'NA': 
#             print l
        for item in l:                
            if item == 'NA':
                continue
            else: 
                for frac in item.split(','):
                    tissue = frac.split(':')[0]
                    assoc_gene_dict[gid].add(tissue)
    for gene in assoc_gene_dict.iterkeys():
        j = assoc_gene_dict[gene]
        if len(j)==1:
            cats.add(gene)
    
    assocs = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    positional = collections.defaultdict(lambda: [])
    fileh = open(args.assoc)
    header = fileh.next().strip().split('\t')
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[header.index('transcript_id')]
        gid = line[header.index('gene_id')]
        chrom = meta_dict[tid]['chrom']
        if chrom == 'chrX': chrom = 'chr26'
        if chrom == 'chrY': chrom = 'chr26'
        start = meta_dict[tid]['start']
        pos_val = (int(chrom[3:]), int(start))
        ct_frac_up = line[header.index('ct_frac_up')]
        cn_frac_up = line[header.index('cn_frac_up')]
        clat_frac_up = line[header.index('clat_frac_up')]
        ct_frac_dn = line[header.index('ct_frac_dn')]
        cn_frac_dn = line[header.index('cn_frac_dn')]
        clat_frac_dn = line[header.index('clat_frac_dn')]
        if gid not in cats:
            positional[gid].append((tid, pos_val))
        else:
            if clat_frac_up != 'NA':
                tissue, val = clat_frac_up.split(':')
                assocs[tissue][gid].append((tid, float(val)))
            elif cn_frac_up != 'NA':         
                tissue, val = cn_frac_up.split(':')
                assocs[tissue][gid].append((tid, float(val)))
            elif ct_frac_up != 'NA':
                tissue, val = ct_frac_up.split(':')
                assocs[tissue][gid].append((tid, float(val)))
            elif clat_frac_dn != 'NA':
                tissue, val = clat_frac_dn.split(':')
                assocs[tissue][gid].append((tid, -1*float(val)))
            elif cn_frac_dn != 'NA':
                tissue, val = cn_frac_dn.split(':')
                assocs[tissue][gid].append((tid, -1*float(val)))
            elif ct_frac_dn != 'NA':
                tissue, val = ct_frac_dn.split(':')
                assocs[tissue][gid].append((tid, -1*float(val)))
            
    
    sorted_lists = collections.defaultdict(lambda: [])
    for tissue in assocs: 
        d = assocs[tissue]
        for gene in d.iterkeys():
            d[gene] = sorted(d[gene], key=itemgetter(1), reverse=True)
        sorted_lists[tissue] = sorted(d.items(), key=lambda e: e[1][0][1], reverse=True)

    suff_counter = collections.defaultdict(lambda: 1)
    renamed = collections.defaultdict(lambda: 'null')
    name_dict = collections.defaultdict(lambda: 'NA')
    type_dict = collections.defaultdict(lambda: 'NA')
    category_dict = collections.defaultdict(lambda: 'NA')
    dir_dict = collections.defaultdict(lambda: 'NA')
    name_list = []
    nonsinglets = []

    for tissue in sorted_lists.iterkeys():
        pre = namer(tissue) + 'AT'
        for g, ts in sorted_lists[tissue]:
            rank = 1
            suff_check = str(rank)
            #make sure name that's generated is not already in the reference list
            name_check = pre + suff_check
            j=0
            while name_check in hugos: 
                j +=1
                suff_check = str(rank + j)
                name_check = pre + suff_check
            for t in ts:
                t_id = t[0]
                meta = meta_dict[t_id]
                category = meta['category']
                ref_name = set(meta['ref_gene_name'].upper().split(','))
                if t[1] >= 0:
                        dir_dict[t_id] = 'up'
                else:
                    dir_dict[t_id] = 'dn'
                if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
                    ref_name = list(ref_name.intersection(lncrna))[0]
                    name = ref_name + '.' + str(suff_counter[ref_name])
                    renamed[t_id] = '0'
                    if suff_counter[ref_name] == 2:
                        nonsinglets.append(ref_name)
                    suff_counter[ref_name] += 1                    
                    name_dict[t_id] = name
                    name_list.append(name)
                else: 
                    suff = suff_check + '.' + str(suff_counter[name_check])
                    name = pre + suff
                    renamed[t_id] = '1'
                    if suff_counter[name_check] == 2:
                        nonsinglets.append(name_check)
                    suff_counter[name_check] +=1
                    hugos.add(name_check)
                    rank +=1
                    name_dict[t_id] = name
                    name_list.append(name)
        
        
    #sort positional dict:
    positionals = sorted(positional.items(), key=lambda e: (e[1][0][1][0], e[1][0][1][1]))
    pre = 'CAT'
    for gene, ts in positionals:
        rank = 1
        suff_check = str(rank)
        name_check = pre + suff_check
        j=0
        while name_check in hugos: 
            j +=1
            suff_check = str(rank + j)
            name_check = pre + suff_check
        for t in ts:
            t_id = t[0]
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
                name_dict[t_id] = name
                name_list.append(name)
            else: 
                suff = suff_check + '.' + str(suff_counter[name_check])
                name = pre + suff
                renamed[t_id] = '1'
                if suff_counter[name_check] == 2:
                    nonsinglets.append(name_check)
                suff_counter[name_check] +=1
                hugos.add(name_check)
                rank +=1
                name_dict[t_id] = name
                name_list.append(name)
    
    #do the hiclincs 
    hiclinc_dict = collections.defaultdict(lambda: [])
    positional_set = set(positional.keys())
    for tid in hiclincs: 
        gid = meta_dict[tid]['gene_id']
        if gid in cats: continue
        if gid in positional_set: continue
        chrom = meta_dict[tid]['chrom']
        if chrom == 'chrX': chrom = 'chr26'
        if chrom == 'chrY': chrom = 'chr26'
        start = meta_dict[tid]['start']
        pos_val = (int(chrom[3:]), int(start))
        hiclinc_dict[gid].append([tid, pos_val])
    
    #sort hiclinc dict:
    hiclinc_list = sorted(hiclinc_dict.items(), key=lambda e: (e[1][0][1][0], e[1][0][1][1]))
    pre = 'HICLINC'
    for gene, ts in hiclinc_list:
        rank = 1
        suff_check = str(rank)
        name_check = pre + suff_check
        j=0
        while name_check in hugos: 
            j +=1
            suff_check = str(rank + j)
            name_check = pre + suff_check
        for t in ts:
            t_id = t[0]
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
                name_dict[t_id] = name
                name_list.append(name)
            else: 
                suff = suff_check + '.' + str(suff_counter[name_check])
                name = pre + suff
                renamed[t_id] = '1'
                if suff_counter[name_check] == 2:
                    nonsinglets.append(name_check)
                suff_counter[name_check] +=1
                hugos.add(name_check)
                rank +=1
                name_dict[t_id] = name
                name_list.append(name)
    
    nonsinglets = set(nonsinglets)
    for key in name_dict.iterkeys():
        name = name_dict[key].split('.')[0]
        if name not in nonsinglets: 
            name_dict[key] = name
    
#     for key, name in name_dict.iteritems():
#         print '\t'.join([key, name])
        
    logging.debug('printing new metadata')
    meta_fh = open(args.meta)
    meta_header = meta_fh.next().strip().split('\t')
    meta_header.append('func_name_final')
    meta_header.append('renamed_final')
    print '\t'.join(meta_header)
    i = 0
    for line in meta_fh:
        i+=1
        if i%10000==0:
            logging.debug('Printed %d lines' % i) 
        line = line.strip().split('\t')
        t_id = line[meta_header.index('transcript_id')]
        func_name = name_dict[t_id]
        if func_name == 'NA':
            continue
        rename = renamed[t_id]
        line.append(func_name)
        line.append(rename )
        print '\t'.join(map(str, line))
            
            
            
            
            
#         
#         
#         for tissue in sorted_lists.iterkeys():
#         pre = namer(tissue) + 'AT'
#         for rank, [g, ts] in enumerate(sorted_lists[tissue]):
#             suff_check = str(rank + 1)
#             #make sure name that's generated is not already in the reference list
#             name_check = pre + suff_check
#             j=0
#     
#     
#     print name_dict
#     
    
    
#     suff_counter = collections.defaultdict(lambda: 1)
#     name_dict = collections.defaultdict(lambda: 'NA')
#     type_dict = collections.defaultdict(lambda: 'NA')
#     category_dict = collections.defaultdict(lambda: 'NA')
#     dir_dict = collections.defaultdict(lambda: 'NA')
#     name_list = []
#     nonsinglets = []
    
#     for tissue in sorted_lists.iterkeys():
#         pre = namer(tissue) + 'AT'
#         for rank, [g, ts] in enumerate(sorted_lists[tissue]):
#             suff_check = str(rank + 1)
#             meta = meta_dict[t_id]
#             category = meta['category']
#             ref_name = set(meta['ref_gene_name'].upper().split(','))
#             if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
#                 ref_name = list(ref_name.intersection(lncrna))[0]
#                 name = ref_name + '.' + str(suff_counter[ref_name])
#                 renamed[t_id] = '0'
#                 if suff_counter[ref_name] == 2:
#                     nonsinglets.append(ref_name)
#                 suff_counter[ref_name] += 1
#             
#             print rank, g, ts
            
            
            
            
            
    
    
#         
#             
#         
#     
#     
#     
#     
#     
#     ref_genes = set(ref_genes)
#     hugos = set(hugos)
#     lncrna = set(lncrna)
#     renamed = collections.defaultdict(lambda: 'null')
#     #function to name the clats, cats, lats
#     def type_namer(dict, label):
#         for type in dict.iterkeys():
#             if type == 'pancancer': 
#                 continue
#             g_dict = collections.defaultdict(lambda: [])
#             pre = namer(type)
#             mid = label.upper()
#             for rank, [t_id, g_id, score] in enumerate(sorted(dict[type], key = lambda x: (-x[2]))):
#                 g_dict[g_id].append([t_id, rank, score])
#             g_rank = []
#             for key in g_dict.iterkeys():
#                 rank = g_dict[key][0][1]
#                 g_rank.append([key, rank])
#             for rank_final, [g_id, rank_rough] in enumerate(sorted(g_rank, key = lambda x: (x[1]))):
#                 trans = g_dict[g_id]
#                 annot_counter = 0
#                 for t_id,r, s in trans:
#                     meta = meta_dict[t_id]
#                     category = meta['category']
#                     ref_name = set(meta['ref_gene_name'].upper().split(','))
#                     if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
#                         annot_counter +=1
#                 lnc_counter = len(trans) - annot_counter
#                 
#                 #make sure name that's generated is not already in the reference list
#                 suff_check = str(rank_final + 1)
#                 name_check = pre + mid + suff_check
#                 j=0
#                 while name_check in hugos: 
#                     j +=1
#                     suff_check = str(rank_final + 1 + j)
#                     name_check = pre + mid + suff_check
#                 if lnc_counter > 0: 
#                     hugos.add(name_check)
#                 
#                 for i,[t_id,r, s] in enumerate(trans): 
#                     #check to see if it is an already named lncrna (from hugo) 
#                     meta = meta_dict[t_id]
#                     category = meta['category']
#                     ref_name = set(meta['ref_gene_name'].upper().split(','))
#                     if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
#                         ref_name = list(ref_name.intersection(lncrna))[0]
#                         name = ref_name + '.' + str(suff_counter[ref_name])
#                         renamed[t_id] = '0'
#                         if suff_counter[ref_name] == 2:
#                             nonsinglets.append(ref_name)
#                         suff_counter[ref_name] += 1
#                     else: 
#                         suff = suff_check + '.' + str(suff_counter[name_check])
#                         if suff_counter[name_check] == 2:
#                             nonsinglets.append(name_check)
#                         suff_counter[name_check] +=1
#                         name = pre + mid + suff
#                         renamed[t_id] = '1'
#                     if s >= 0:
#                         dir_dict[t_id] = 'up'
#                     else:
#                         dir_dict[t_id] = 'dn'
#                     name_dict[t_id] = name
#                     type_dict[t_id] = type
#                     if t_id in clats:
#                         category_dict[t_id] = 'clat'
#                     else: 
#                         category_dict[t_id] = label.lower()
#                     name_list.append(name)
#                     
#     
#    
#     #name the clats, cats, lats
#     type_namer(clat_dict, 'cat')
#     type_namer(cat_dict, 'cat')
#     type_namer(lat_dict, 'at')
#         
#         
#     #name the hiclincs
#     g_dict = collections.defaultdict(lambda: [])
#     mid = 'HICLINC'
#     for rank, [t_id, g_id, chrom, start] in enumerate(sorted(hiclincs, key = lambda x: (x[2], x[3]))):
#         g_dict[g_id].append([t_id, rank])
#     g_rank = []
#     for key in g_dict.iterkeys():
#         rank = g_dict[key][0][1]
#         g_rank.append([key, rank])
#     for rank_final, [g_id, rank_rough] in enumerate(sorted(g_rank, key = lambda x: (x[1]))):
#         trans = g_dict[g_id]
#         annot_counter = 0
#         for t_id,r in trans:
#             meta = meta_dict[t_id]
#             category = meta['category']
#             ref_name = set(meta['ref_gene_name'].upper().split(','))
#             if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
#                 annot_counter +=1
#         lnc_counter = len(trans) - annot_counter
#         
#         #make sure name that's generated is not already in the reference list
#         suff_check = str(rank_final + 1)
#         name_check = mid + suff_check
#         j=0
#         while name_check in hugos: 
#             j +=1
#             suff_check = str(rank_final + 1 + j)
#             name_check = mid + suff_check
#         if lnc_counter > 0: 
#             hugos.add(name_check)
#         
#         for i,[t_id,r] in enumerate(trans): 
#             if t_id == 'T135389':
#                 continue
#             #check to see if it is an already named lncrna (from hugo) 
#             meta = meta_dict[t_id]
#             category = meta['category']
#             ref_name = set(meta['ref_gene_name'].upper().split(','))
#             if category in ['same_strand', 'read_through'] and len(ref_name.intersection(lncrna))>0:
#                 ref_name = list(ref_name.intersection(lncrna))[0]
#                 name = ref_name + '.' + str(suff_counter[ref_name])
#                 renamed[t_id] = '0'
#                 if suff_counter[ref_name] == 2:
#                     nonsinglets.append(ref_name)
#                 suff_counter[ref_name] += 1
#             else: 
#                 suff = suff_check + '.' + str(suff_counter[name_check])
#                 if suff_counter[name_check] == 2:
#                     nonsinglets.append(name_check)
#                 suff_counter[name_check] +=1
#                 name = mid + suff
#                 renamed[t_id] = '1'
#             name_dict[t_id] = name
#             type_dict[t_id] = 'hiclinc'
#             category_dict[t_id] = 'hiclinc'
#             name_list.append(name)
# 
#     nonsinglets = set(nonsinglets)
#     for key in name_dict.iterkeys():
#         name = name_dict[key].split('.')[0]
#         if name not in nonsinglets: 
#             name_dict[key] = name
#     
#         
#     logging.debug('printing new metadata')
#     meta_fh = open(args.meta)
#     meta_header = meta_fh.next().strip().split('\t')
#     meta_header.append('func_name')
#     meta_header.append('func_type')
#     meta_header.append('func_cat')
#     meta_header.append('func_dir')
#     meta_header.append('renamed')
#     print '\t'.join(meta_header)
#     i = 0
#     ref_genes = []
#     for line in meta_fh:
#         i+=1
#         if i%10000==0:
#             logging.debug('Printed %d lines' % i) 
#         line = line.strip().split('\t')
#         t_id = line[meta_header.index('transcript_id')]
#         func_name = name_dict[t_id]
#         func_type = type_dict[t_id]
#         func_cat = category_dict[t_id]
#         func_dir = dir_dict[t_id]
#         rename = renamed[t_id]
#         line.append(func_name)
#         line.append(func_type)
#         line.append(func_cat)
#         line.append(func_dir)
#         line.append(rename )
#         print '\t'.join(map(str, line))
    return 1
        
        
        
if __name__ == '__main__': 
    sys.exit(main())
