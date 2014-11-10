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
    ct_dict_up = collections.defaultdict(lambda: set())
    ct_dict_dn = collections.defaultdict(lambda: set())
    cn_dict_up = collections.defaultdict(lambda: set())
    cn_dict_dn = collections.defaultdict(lambda: set())
    clat_dict_up = collections.defaultdict(lambda: set())
    clat_dict_dn = collections.defaultdict(lambda: set())
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
    stats_dict = collections.defaultdict(lambda: 0)
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
        key = (analysis, can_type, dir, tid)
        stats_dict[key] = frac
        
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
                if dir == 'up': 
                    ct_dict_up[tid].add(can_type)
                else: 
                    ct_dict_dn[tid].add(can_type)
                ct_trans.add(tid)
                ct_gene.add(gid)
                tissues.append(can_type)
                
        if analysis == 'cn': 
            if frac > CN_FRAC_CUTOFF and fdr < FDR_CUTOFF:
                if dir == 'up':
                    cn_dict_up[tid].add(can_type)
                else: 
                    cn_dict_dn[tid].add(can_type)
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
            clat_dict_up[tid] = cn & ct
            clat_trans.add(tid)
            clat_gene.add(gid)
            tissues = tissues + list(cn&ct)
    for tid in clat_cand_dict_dn.iterkeys():
        cn = clat_cand_dict_dn[tid]
        ct = clat_ct_cand_dict_dn[tid]
        gid = gene_dict[tid]
        if len(cn & ct) > 0:
            clat_dict_dn[tid] = cn & ct
            clat_trans.add(tid)
            clat_gene.add(gid)
            tissues = tissues + list(cn&ct)
    
    
    print '\t'.join(['transcript_id', 'gene_id', 
                     'ct_up', 'ct_frac_up',
                     'cn_up', 'cn_frac_up',
                     'clat_up', 'clat_frac_up',
                     'ct_dn', 'ct_frac_dn', 
                     'cn_dn', 'cn_frac_dn',
                     'clat_dn', 'clat_frac_dn'])
    for tid in ct_trans|cn_trans|clat_trans:
        gid = gene_dict[tid]
        ct_up_tissue_types = ','.join(ct_dict_up[tid])
        ct_dn_tissue_types = ','.join(ct_dict_dn[tid])
        cn_up_tissue_types = ','.join(cn_dict_up[tid])
        cn_dn_tissue_types = ','.join(cn_dict_dn[tid])
        clat_up_tissue_types = ','.join(clat_dict_up[tid])
        clat_dn_tissue_types = ','.join(clat_dict_dn[tid])
        
        if ct_up_tissue_types == '':
            ct_up_tissue_types = 'NA'
            ct_up_frac = ['NA']
        else: 
            ct_up_frac = []
            for tissue in ct_up_tissue_types.split(','):
                key = ('ct', tissue, 'up', tid)
                val = stats_dict[key]
                ct_up_frac.append('%s:%s' % (tissue, val))
                
        if cn_up_tissue_types == '':
            cn_up_tissue_types = 'NA'
            cn_up_frac = ['NA']
        else: 
            cn_up_frac = []
            for tissue in cn_up_tissue_types.split(','):
                key = ('cn', tissue, 'up', tid)
                val = stats_dict[key]
                cn_up_frac.append('%s:%s' % (tissue, val))
        
        if clat_up_tissue_types == '':
            clat_up_tissue_types = 'NA'
            clat_up_frac = ['NA']
        else:
            clat_up_frac = [] 
            for tissue in clat_up_tissue_types.split(','):
                key = ('ct', tissue, 'up', tid)
                ct_frac = stats_dict[key]
                key = ('cn', tissue, 'up', tid)
                cn_frac = stats_dict[key]
                val = (ct_frac + cn_frac)/2.0
                clat_up_frac.append('%s:%s' % (tissue, val))
        
        if ct_dn_tissue_types == '':
            ct_dn_tissue_types = 'NA'
            ct_dn_frac = ['NA']
        else: 
            ct_dn_frac = []
            for tissue in ct_dn_tissue_types.split(','):
                key = ('ct', tissue, 'dn', tid)
                val = stats_dict[key]
                ct_dn_frac.append('%s:%s' % (tissue, val))
                
        if cn_dn_tissue_types == '':
            cn_dn_tissue_types = 'NA'
            cn_dn_frac = ['NA']
        else: 
            cn_dn_frac = []
            for tissue in cn_dn_tissue_types.split(','):
                key = ('cn', tissue, 'dn', tid)
                val = stats_dict[key]
                cn_dn_frac.append('%s:%s' % (tissue, val))
        
        if clat_dn_tissue_types == '':
            clat_dn_tissue_types = 'NA'
            clat_dn_frac = ['NA']
        else:
            clat_dn_frac = [] 
            for tissue in clat_dn_tissue_types.split(','):
                key = ('ct', tissue, 'dn', tid)
                ct_frac = stats_dict[key]
                key = ('cn', tissue, 'dn', tid)
                cn_frac = stats_dict[key]
                val = (ct_frac + cn_frac)/2.0
                clat_dn_frac.append('%s:%s' % (tissue, val))
        
        
        lineo = [tid, gid, 
                 ct_up_tissue_types, ','.join(ct_up_frac),
                 cn_up_tissue_types, ','.join(cn_up_frac), 
                 clat_up_tissue_types, ','.join(clat_up_frac),
                 ct_dn_tissue_types, ','.join(ct_dn_frac), 
                 cn_dn_tissue_types, ','.join(cn_dn_frac), 
                 clat_dn_tissue_types, ','.join(clat_dn_frac)]
#         if len(clat_dict[tid]) >= 1 and len(cn_dict[tid]) >=0:
#             print tid
#             print clat_dict[tid]
#             print cn_dict[tid]
#             print 'tit'
#             cn_tally+=1
        print '\t'.join(map(str,lineo))
#     print cn_tally
    
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
