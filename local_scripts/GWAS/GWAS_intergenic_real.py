'''
Created on Jan 16, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import numpy as np
import collections
import subprocess
import multiprocessing
import shutil

from assemblyline.lib.base import which
from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.gtf import sort_gtf
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.transcript import parse_gtf, strand_int_to_str


'''
1) will intersect real assembly results with gwas snps to check for overlap
2) prints info for that intersection
''' 

intergenic_assembly_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.intergenic.sorted.bed'
snp_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed'
excl_file = '/mctp/projects/mitranscriptome/gwas/bed_files/ref.merged.gaps.sorted.bed'
chrom_sizes_file = '/mctp/projects/rnaseq/version_001_2013-01-14/references/Oncoseq_hg19_Linux_x86-64/Annotation/ChromInfo.txt'
gtf_file = '/mctp/projects/mitranscriptome/gwas/gtfs/intergenic.gtf'
gwas_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/rand_snps2.sorted.bed'

RSIDCOL = 3
x = 0



class Interval(object):
    def __init__(self):
        self.gene_id = None
        self.chrom = None
        self.start = None
        self.end = None

def get_gene_intervals(transcripts):
    gene_map = collections.defaultdict(lambda: Interval())
    for t in transcripts:
        gene_id = t.attrs["gene_id"]
        if gene_id not in gene_map:
            g = Interval()
            g.gene_id = gene_id
            g.chrom = t.chrom
            g.start = t.start
            g.end = t.end
            gene_map[gene_id] = g
        else:
            g = gene_map[gene_id]
        # update interval
        g.start = min(g.start, t.start)
        g.end = max(g.end, t.end)
    for g in gene_map.itervalues():
        yield g

def write_bed(chrom, name, strand, score, exons, flank, chrom_length):  
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    chr_len = chrom_length[chrom]
    tx_start = exons[0].start
    tx_start = max(0, (tx_start - flank))
    tx_end = exons[-1].end    
    tx_end = min(chr_len, (tx_end + flank))
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)        
    # make bed fields
    fields = [chrom, 
              str(tx_start), 
              str(tx_end),
              str(name),
              str(score),
              strand_int_to_str(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str,block_sizes)) + ',',
              ','.join(map(str,block_starts)) + ',']
    return fields


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", dest = 'assembly_bed',
                    default = intergenic_assembly_bed,
                    help = 'Assembly file used for shuffling and snp overlap intersection')
    parser.add_argument("--snps", dest = 'snps',
                    default = snp_bed,
                    help = 'SNP universe bed file')
    parser.add_argument("--excl", dest = 'excl',
                    default = excl_file,
                    help = 'Exclusion file used for shuffling')
    parser.add_argument("--chrom", dest = 'chrom',
                    default = chrom_sizes_file,
                    help = 'Chrom size file used for shuffling')
    parser.add_argument("--gtf", dest = 'gtf',
                    default = gtf_file,
                    help = 'GTF file used to generate shuffle (should match assembly_bed)')
    parser.add_argument("--gwas", dest = 'gwas',
                    default = gwas_bed,
                    help = 'GWAS bed file file used for intersection')
    parser.add_argument("--flank", dest = 'flank',
                    default = 0,
                    help = 'number of flanking bases to add to bed files')
    args = parser.parse_args()
    
    args.flank = int(args.flank)
    
    logging.info('Output is printed to stdout, to save use \'>\' <filename>')
    
    # check command line parameters
    if which('bedtools') is None:
        parser.error('bedtools binary not found in PATH')
    if not os.path.exists(chrom_sizes_file):
        parser.error('chrom sizes file %s not found' % (chrom_sizes_file))

    if not os.path.isdir('GWAS_TMPS'):
        os.mkdir('GWAS_TMPS')
    
    prefix = 'GWAS_TMPS'    
    gene_intervals_file = os.path.join(prefix, 'gene_intervals.bed')
    intersect_file = os.path.join(prefix, 'intersect.txt')
    assembly_flank = os.path.join(prefix, 'flank.bed')
    
    output_file = 'gwas_intergenic_null.txt'
    
        
    logging.info('Parsing GTF file')
    with open(gene_intervals_file, 'w') as f:
        for locus_transcripts in parse_gtf(open(gtf_file)):
            # find borders of locus
            locus_chrom = locus_transcripts[0].chrom
            locus_start = min(t.start for t in locus_transcripts)
            locus_end = max(t.end for t in locus_transcripts)
#             logging.debug("[LOCUS] %s:%d-%d %d transcripts" % 
#                           (locus_chrom, locus_start, locus_end, 
#                            len(locus_transcripts)))
            for g in get_gene_intervals(locus_transcripts):
                print >>f, '\t'.join(map(str, [g.chrom, g.start, g.end, g.gene_id]))   
                
    #apply flank to the bed file 
    #read chrom file to make sure flanks added do not enter chrom ends
    chrom_length = {}
    for line in open(chrom_sizes_file): 
        line = line.strip().split('\t')
        chr = line[0]
        length = line[1]
        chrom_length[chr] = length
    with open(assembly_flank, 'w') as f:
        for line in open(args.assembly_bed):
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            chr_len = chrom_length[chr]
            start = max(0, (start - args.flank))
            end = min(chr_len, (end + args.flank))
            line[1] = start
            line[2] = end
            print >> f, '\t'.join(map(str, line))
    
    
    #GWAS snps
    #do intersections for real data and report number of overlapping GWAS snps 
    logging.info('Intersecting assembly with GWAS snps')
    args_int = ['bedtools', 'intersect', 
            '-a', args.gwas,
            '-b', assembly_flank,
            '-wa',
            '-wb']
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args_int, stdout=fileh)
    #count number of SNPs caught
    snps = set()
    
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    gwas_overlap = len(snps)

    #snp universe
    #do intersections for real data and report number of overlapping snps in snp universe 
    logging.info('Intersecting assembly with snp universe')
    args_int = ['bedtools', 'intersect', 
            '-a', args.snps,
            '-b', assembly_flank,
            '-wa',
            '-wb',
            '-sorted']
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args_int, stdout=fileh)
    #count number of SNPs caught
    snps = set()
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    snp_overlap = len(snps)
    frac_real = float(gwas_overlap)/snp_overlap
    logging.info('%d GWAS snps overlap compendia genes'  % gwas_overlap)
    logging.info('%d snps (from \"snp universe\") overlap compendia genes' % snp_overlap)
    logging.info('Frac: %f' % frac_real)
    
    
    print '\t'.join(map(str, [args.flank, gwas_overlap, snp_overlap, frac_real]))
    

    
#     valso.close()
    shutil.rmtree(prefix)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())



