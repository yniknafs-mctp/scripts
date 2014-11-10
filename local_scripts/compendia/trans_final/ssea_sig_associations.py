'''
Created on Aug 15, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import collections

'''
identify all of the significant cancer associations for the lncRNAs 
takes a ssea query file and the metadata
'''

PREC_CUTOFF = .5
CT_FRAC_CUTOFF = .975
CT_NOPREC_CUTOFF = .995
CN_FRAC_CUTOFF = .978285
CLAT_CUTOFF = .95
FDR_CUTOFF = .00001

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

def translate(set_name):
    fields = set_name.split('_')
    prog = fields[0]
    if fields[1] == 'versus':
        analysis = 'cn'
        can_type = '_'.join(fields[3:])
    else: 
        analysis = 'ct'
        can_type = '_'.join(fields[2:])
    
    if can_type == 'skeletal_muscle':
        prog = 'cancer'
    if can_type =='heart':
        prog = 'cancer'
    if can_type == 'embryonic_stem_cells':
        prog = 'cancer'
    if can_type == 'breast_carcinoma':
        can_type = 'breast'
    if can_type == 'cervical_carcinoma':
        can_type = 'cervical'
    if can_type == 'colorectal_carcinoma':
        can_type = 'colorectal'
    if can_type == 'glioblastoma_multiforme_gbm':
        can_type = 'gbm'
    if can_type == 'head_and_neck_carcinoma':
        can_type = 'head_neck'
    if can_type == 'hepatocellular_carcinoma':
        can_type = 'liver'
    if can_type == 'lower_grade_glioma_lgg':
        can_type = 'lgg'
    if can_type == 'lung_adenocarcinoma':
        can_type = 'luad'
    if can_type == 'lung_squamous_cell_carcinoma':
        can_type = 'lusc'
    if can_type == 'pancreatic_carcinoma':
        can_type = 'pancreatic'
    if can_type == 'prostate_carcinoma':
        can_type = 'prostate'
    if can_type == 'renal_cell_carcinoma_chromophobe':
        can_type = 'kich'
    if can_type == 'renal_clear_cell_carcinoma':
        can_type = 'kirc'
    if can_type == 'renal_papillary_cell_carcinoma':
        can_type = 'kirp'
    if can_type == 'stomach_adenocarcinoma':
        can_type = 'stomach'
    if can_type == 'uterine_endometrial_carcinoma':
        can_type = 'uterine'
    if can_type == 'head_and_neck':
        can_type = 'head_neck'
    if can_type == 'kidney_chromophobe':
        can_type = 'kich'
    if can_type == 'kidney_clear_cell':
        can_type = 'kirc'
    if can_type == 'kidney_papillary_cell':
        can_type = 'kirp'
    if can_type == 'lung_adenocarcinoma':
        can_type = 'luad'
    if can_type == 'lung_squamous':
        can_type = 'lusc'

    return prog, analysis, can_type
    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("ssea")
    parser.add_argument("metadata")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Reading metadata')
    fileh = open(args.metadata)
    header = fileh.next().strip().split('\t')
    gene_dict = {}
    tcat_dict = {}
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[header.index("transcript_id")]
        gid = line[header.index("gene_id")]
        tcat = line[header.index("tcat")]
        gene_dict[tid] = gid
        tcat_dict[tid] = tcat
    
    
    logging.info('Reading SSEA data')
    fileh = open(args.ssea)
    header = fileh.next().strip().split("\t")
    ct_dict = collections.defaultdict(lambda: set())
    cn_dict = collections.defaultdict(lambda: set())
    clat_dict = collections.defaultdict(lambda: set())
    clat_cand_dict_up = collections.defaultdict(lambda: set())
    clat_ct_cand_dict_up = collections.defaultdict(lambda: set()) 
    clat_cand_dict_dn = collections.defaultdict(lambda: set())
    clat_ct_cand_dict_dn = collections.defaultdict(lambda: set()) 
    ct_trans = set()
    cn_trans = set()
    clat_trans = set()
    ct_gene = set()
    cn_gene = set()
    clat_gene = set()
    tissues = []
    for line in fileh: 
        line = line.strip().split("\t")
        tid = line[header.index("transcript_id")]
        gid = gene_dict[tid]
        tcat = tcat_dict[tid]
        if tcat not in ['lncrna', 'tucp']:
            continue
        set_name = line[header.index('ss_compname')]
        frac = float(line[header.index('frac')])
        if frac > 0: 
            dir = 'up'
        else: 
            dir = 'dn'
        frac = abs(frac)
        prec = float(line[header.index('prec')])
        fdr = float(line[header.index('fdr')])
        prog, analysis, can_type = translate(set_name)
        if prog == 'normal':
            continue
        if analysis == 'ct':
            if can_type == 'all':
                continue 
            if can_type in no_lat:
                    if frac > CLAT_CUTOFF and prec > PREC_CUTOFF:
                        if dir == 'up':
                            clat_ct_cand_dict_up[tid].add(can_type)
                        else: 
                            clat_ct_cand_dict_dn[tid].add(can_type)
                    continue
            if (frac > CT_FRAC_CUTOFF and prec > PREC_CUTOFF) or frac > CT_NOPREC_CUTOFF:
                ct_dict[tid].add(can_type)
                ct_trans.add(tid)
                ct_gene.add(gid)
                tissues.append(can_type)
                
        if analysis == 'cn': 
            if frac > CN_FRAC_CUTOFF and fdr < FDR_CUTOFF:
                cn_dict[tid].add(can_type)
                cn_trans.add(tid)
                cn_gene.add(gid)
                tissues.append(can_type)
            if frac > CLAT_CUTOFF and fdr < FDR_CUTOFF:
                if dir == 'up':
                    clat_cand_dict_up[tid].add(can_type)
                else: 
                    clat_cand_dict_dn[tid].add(can_type)
                
        
    for tid in clat_cand_dict_up.iterkeys():
        cn = clat_cand_dict_up[tid]
        ct = clat_ct_cand_dict_up[tid]
        gid = gene_dict[tid]
        if len(cn & ct) > 0:
            clat_dict[tid] = cn & ct
            clat_trans.add(tid)
            clat_gene.add(gid)
            tissues = tissues + list(cn&ct)
    for tid in clat_cand_dict_dn.iterkeys():
        cn = clat_cand_dict_dn[tid]
        ct = clat_ct_cand_dict_dn[tid]
        gid = gene_dict[tid]
        if len(cn & ct) > 0:
            clat_dict[tid] = cn & ct
            clat_trans.add(tid)
            clat_gene.add(gid)
            tissues = tissues + list(cn&ct)
    
    
#     print '\t'.join(['transcript_id', 'ct', 'cn', 'clat'])
    cn_tally = 0
    for tid in ct_trans|cn_trans|clat_trans:
        ct_tissue_types = ','.join(ct_dict[tid])
        cn_tissue_types = ','.join(cn_dict[tid])
        clat_tissue_types = ','.join(clat_dict[tid])
        if ct_tissue_types == '':
            ct_tissue_types = 'NA'
        if cn_tissue_types == '':
            cn_tissue_types = 'NA'
        if clat_tissue_types == '':
            clat_tissue_types = 'NA'
        lineo = [tid, ct_tissue_types, cn_tissue_types, clat_tissue_types]
        if len(clat_dict[tid]) >= 1 and len(cn_dict[tid]) >=0:
            print tid
            print clat_dict[tid]
            print cn_dict[tid]
            print 'tit'
            cn_tally+=1
#         print '\t'.join(lineo)
    print cn_tally
    
#     for gid in ct_trans|cn_trans|clat_trans:
#         print gid
    
        
#     logging.debug('CT: %d' % len(ct_trans))
#     logging.debug('CN: %d' % len(cn_trans))
#     logging.debug('CLAT: %d' % len(clat_trans))
#     logging.debug('All: %d' % len(clat_trans|ct_trans|cn_trans))
#     logging.debug('All_g: %d' % len(clat_gene|ct_gene|cn_gene))
#     logging.debug(collections.Counter(tissues))
#     logging.debug(len(collections.Counter(tissues)))
#     logging.debug('G037188' in ct_trans)
#     logging.debug('G037188' in cn_trans)
#     logging.debug('G037188' in clat_trans)
#     
        
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
