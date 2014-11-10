'''
Created on Aug 5, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import collections




PROSTATE = [
            'prostate_lncap_si_3812',
            'prostate_lncap_si_3811'
#             'prostate_du145_si_3193_ga2'
            
            ]

BREAST = [
            'idea_breast_mcf7_50bp_pe',
            'idea_breast_mcf7_100bp_sr',
            'encode_breast_mcf7_caltech_rep3',
            'encode_breast_mcf7_caltech_rep2',
            'encode_breast_mcf7_caltech_rep1'
#             'idea_breast_mda_mb_231_50bp_pe',
#             'idea_breast_mda_mb_231_100bp_sr',
#             'breast_mda_mb_231_si_2870'
          ]

LUNG = [
            'lung_cl_a549_si_4870',
            'encode_lung_cl_a549_cshl_rep2',
            'encode_lung_cl_a549_cshl_rep1'
#             'lung_cl_nci_h838_si_4622',
        ]

CELL_LINE_EXPR_CO = 1
SSEA_CUTOFF = .97
SSEA_COMBO_CUTOFF = .95

def ssea_parse(ssea_dict, tid):
    sig_cats = []
    ct_b = ssea_dict['cancer_type_breast_carcinoma'][tid]
    cn_b = ssea_dict['cancer_versus_normal_breast'][tid]
    b_sign = ct_b*cn_b
    ct_b = abs(ct_b)
    cn_b = abs(cn_b)
    ct_p = ssea_dict['cancer_type_prostate_carcinoma'][tid]
    cn_p = ssea_dict['cancer_versus_normal_prostate'][tid]
    p_sign = ct_p*cn_p
    ct_p = abs(ct_p)
    cn_p = abs(cn_p)
    ct_l = ssea_dict['cancer_type_lung_adenocarcinoma'][tid]
    cn_l = ssea_dict['cancer_versus_normal_lung_adenocarcinoma'][tid]
    l_sign = ct_l*cn_l
    ct_l = abs(ct_l)
    cn_l = abs(cn_l)
    if ct_b > SSEA_COMBO_CUTOFF and cn_b> SSEA_COMBO_CUTOFF:
        if b_sign > 0:
            sig_cats.append('breast')
    if ct_p > SSEA_COMBO_CUTOFF and cn_p> SSEA_COMBO_CUTOFF:
        if p_sign > 0:
            sig_cats.append('prostate')
    if ct_l > SSEA_COMBO_CUTOFF and cn_l> SSEA_COMBO_CUTOFF:
        if l_sign > 0:
            sig_cats.append('luad')
    if cn_b > SSEA_CUTOFF:
        sig_cats.append('breast')
    if cn_p > SSEA_CUTOFF:
        sig_cats.append('prostate')
    if cn_l > SSEA_CUTOFF:
        sig_cats.append('luad')
    return set(sig_cats)

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata")
    parser.add_argument("cell_line_expr")
    parser.add_argument("ssea")
    parser.add_argument("sequences")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Reading cell line expr data')
    cl_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    fileh = open(args.cell_line_expr)
    header = fileh.next().strip().split('\t')
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0]
        for lib in header[1:]:
            cl_dict[lib][tid] = float(line[header.index(lib)])
    
    logging.info('Reading ssea results data')
    ssea_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    fileh = open(args.ssea)
    header = fileh.next().strip().split('\t')
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0]
        set_name = line[header.index('ss_compname')] 
        frac = line[header.index('frac')]
        ssea_dict[set_name][tid] = float(frac)
    
    logging.info('Reading sequence data')
    seq_dict = collections.defaultdict(lambda: 'na')
    fileh = open(args.sequences)
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[0].split('|')[1]
        seq = line[2].upper().replace('|','[]')
        seq_dict[tid] = seq
    
    fileh = open(args.metadata)
    header = fileh.next().strip().split('\t')
    monoexonic = []
    exon_dict = collections.defaultdict(lambda: [])
    for line in fileh: 
        line = line.strip().split('\t')
        tid = line[header.index('transcript_id')]
        gid = line[header.index('gene_id')]
        num_exons = line[header.index('num_exons')]
        exon_dict[gid].append(num_exons)
    for key in exon_dict.iterkeys():
        exons = exon_dict[key]
        if len(set(exons)) == 1 and list(set(exons))[0] == '1':
            monoexonic.append(key)
    monoexonic = set(monoexonic)
      
    fileh = open(args.metadata)
    header = fileh.next().strip().split('\t')
    headero = [
               'transcript_id',
               'gene_id',
               'func_name',
               'chrom',
               'start',
               'end',
               'num_exon',
               'monoexonic_gene',
               'tissues'
               ]
    for cl in cl_dict.iterkeys():
        if cl in BREAST or cl in PROSTATE or cl in LUNG:
            headero.append(cl)
            logging.debug(cl)
    headero.append('sequence')
    print '\t'.join(headero)
    
    t_counter_b = 0
    t_counter_p = 0
    t_counter_l = 0
    t_counter_mono_b = 0
    t_counter_mono_p = 0
    t_counter_mono_l = 0
    g_list_b = []
    g_list_p = []
    g_list_l = []
    g_list = []
    g_list_mono_b = []
    g_list_mono_p = []
    g_list_mono_l = []
    g_list_mono = []
    bl = []
    bp = []
    pl = []
    bpl = []
    na = 0
    mono_gene_counter = []
    
    for line in fileh:
        line = line.strip().split('\t')
        tid = line[header.index('transcript_id')]
        gid = line[header.index('gene_id')]
        chrom = line[header.index('chrom')]
        start = line[header.index('start')]
        end = line[header.index('end')]
        num_exon = line[header.index('num_exons')]
        tcat = line[header.index('tcat')]
        func_name = line[header.index('func_name')]
        seq = seq_dict[tid]
#         line = line.strip().split('\t')
#         tid = line[0].split('|')[1]
#         seq = line[2].upper().replace('|','[]')
#         seq_dict[tid] = seq
            

        if tcat not in ['lncrna', 'tucp']:
            continue
        func_type = line[header.index('func_type')]
        sig_cats = ssea_parse(ssea_dict, tid)

        if 'breast' in sig_cats:
            cl_test = 0
            for cl in BREAST:
                if cl_dict[cl][tid] > CELL_LINE_EXPR_CO:
                    cl_test +=1
        elif 'prostate' in sig_cats:
            cl_test = 0
            for cl in PROSTATE:
                if cl_dict[cl][tid] > CELL_LINE_EXPR_CO:
                    cl_test +=1
        elif 'luad' in sig_cats:
            cl_test = 0
            for cl in LUNG:
                if cl_dict[cl][tid] > CELL_LINE_EXPR_CO:
                    cl_test +=1
        else: 
            continue
        if cl_test > 0:
#             if func_type == 'NA': 
#                 na+=1
#                 logging.debug('NA: %s; %s; exons: %s' % (tid, str(sig_cats), num_exon))
#                 continue
            if 'breast' in sig_cats:
                t_counter_b +=1
                g_list_b.append(gid)
                g_list.append(gid)
                if num_exon == '1':
                    t_counter_mono_b +=1
                    g_list_mono_b.append(gid)
                    g_list_mono.append(gid)
            if 'prostate' in sig_cats:
                t_counter_p +=1
                g_list_p.append(gid)
                g_list.append(gid)
                if num_exon == '1':
                    t_counter_mono_p +=1
                    g_list_mono_p.append(gid)
                    g_list_mono.append(gid)
            if 'luad' in sig_cats:
                t_counter_l +=1
                g_list_l.append(gid)
                g_list.append(gid)
                if num_exon == '1':
                    t_counter_mono_l +=1
                    g_list_mono_l.append(gid)
                    g_list_mono.append(gid)
            if 'luad' in sig_cats and 'prostate' in sig_cats:
                pl.append(gid)
            if 'luad' in sig_cats and 'breast' in sig_cats:
                bl.append(gid)
            if 'prostate' in sig_cats and 'breast' in sig_cats:
                bp.append(gid)
            if 'prostate' in sig_cats and 'breast' in sig_cats and 'luad' in sig_cats:
                bpl.append(gid)
            
            if gid in monoexonic: 
                mono_gene = '1'
                mono_gene_counter.append(gid)
            else:
                mono_gene = '0'
            lineo = [
                     tid, 
                     gid,
                     func_name,
                     chrom,
                     start,
                     end,
                     num_exon,
                     mono_gene,
                     ','.join(sig_cats)
                     ]
            for cl in cl_dict.iterkeys():
                if cl in BREAST or cl in PROSTATE or cl in LUNG:
                    lineo.append(cl_dict[cl][tid])
            lineo.append(seq)
            print '\t'.join(map(str,lineo))
    
                
    
    
    logging.info('Breast: %d/%d (%d/%d monoexonic)' % 
                 (t_counter_b, len(set(g_list_b)),
                  t_counter_mono_b, len(set(g_list_mono_b))))
    logging.info('Prostate: %d/%d (%d/%d monoexonic)' % 
                 (t_counter_p, len(set(g_list_p)),
                  t_counter_mono_p, len(set(g_list_mono_p))))
    logging.info('LUAD: %d/%d (%d/%d monoexonic)' % 
                 (t_counter_l, len(set(g_list_l)),
                  t_counter_mono_l, len(set(g_list_mono_l))))
    logging.info('All: %d (%d monoexonic)' % 
                 (len(set(g_list)),len(set(g_list_mono))))
    logging.info('LUAD & breast: %d' % len(set(bl)))
    logging.info('LUAD & prostate: %d' % len(set(pl)))
    logging.info('prostate & breast: %d' % len(set(bp)))
    logging.info('LUAD & breast & prostate: %d' % len(set(bpl)))
    logging.info('NAs: %d' % na)
    logging.info('monoexonic genes: %d' % len(set(mono_gene_counter)&monoexonic))
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
