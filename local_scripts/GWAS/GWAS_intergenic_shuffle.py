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

from assemblyline.lib.base import which
from assemblyline.lib.transcript import parse_gtf
from assemblyline.lib.gtf import sort_gtf
from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.transcript import parse_gtf, strand_int_to_str


'''
1) will shuffle the intergenic genes inside the intergenic space (compared to ref merged bed)
2) will intersect shuffle with gwas snps to check for overlap
3) for each shuffle report number of GWAS snps that are caught
''' 
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

def write_bed(chrom, name, strand, score, exons):
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end    
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

        
    prefix = 'RAND_SNPS_SHUFF'
    intergenic_assembly_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/assembly.intergenic.sorted.bed'
    snp_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/snp_arrays.sorted.bed'
    excl_file = '/mctp/projects/mitranscriptome/gwas/bed_files/ref.merged.gaps.sorted.bed'
    chrom_sizes_file = '/mctp/projects/rnaseq/version_001_2013-01-14/references/Oncoseq_hg19_Linux_x86-64/Annotation/ChromInfo.txt'
    gtf_file = '/mctp/projects/mitranscriptome/gwas/gtfs/intergenic.gtf'
    gwas_bed = '/mctp/projects/mitranscriptome/gwas/bed_files/rand_snps2.sorted.bed'
    # check command line parameters
    if which('bedtools') is None:
        parser.error('bedtools binary not found in PATH')
    if not os.path.exists(chrom_sizes_file):
        parser.error('chrom sizes file %s not found' % (chrom_sizes_file))
    gene_intervals_file = prefix + '.gene_intervals.bed'
    gene_intervals_shuffled_file = prefix + '.gene_intervals.shuffle.bed'
    shuffled_bed_file = prefix + '.shuffle.bed'
    intersect_file = prefix + '.intersect.txt'
    output_file = prefix + '.intergenic_null.txt'
    #keep track of temp files made to delete after finished
    tmp_files = [gene_intervals_file, gene_intervals_shuffled_file, shuffled_bed_file, intersect_file]
    
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
                
    #GWAS snps
    #do intersections for real data and report number of overlapping GWAS snps 
    logging.info('Intersecting assembly with GWAS snps')
    args = ['bedtools', 'intersect', 
            '-a', gwas_bed,
            '-b', intergenic_assembly_bed,
            '-wa',
            '-wb']
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args, stdout=fileh)
    #count number of SNPs caught
    snps = set()
    RSIDCOL = 3
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    gwas_overlap = len(snps)

    #snp universe
    #do intersections for real data and report number of overlapping snps in snp universe 
    logging.info('Intersecting assembly with snp universe')
    args = ['bedtools', 'intersect', 
            '-a', snp_bed,
            '-b', intergenic_assembly_bed,
            '-wa',
            '-wb',
            '-sorted']
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args, stdout=fileh)
    #count number of SNPs caught
    snps = set()
    RSIDCOL = 3
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    snp_overlap = len(snps)
    frac_real = float(gwas_overlap)/snp_overlap
    logging.info('%d GWAS snps overlap compendia genes'  % gwas_overlap)
    logging.info('%d snps (from \"snp universe\") overlap compendia genes' % snp_overlap)
    logging.info('Frac: %f' % frac_real)
    
    
    #loop the shuffle to generate a distribution of nulls for number of snps hit by random intergenic genes
    NUM_SHUFFS = 100
    gwas_vals = []
    universe_vals = []
    valso = open(output_file,'w')
    for x in xrange(NUM_SHUFFS):
#         logging.info('Performing shuffle %d/%d' % (x+1, NUM_SHUFFS))
        args = ['bedtools', 'shuffle', 
                '-excl', excl_file,
                '-i', gene_intervals_file, 
                '-g', chrom_sizes_file]
        with open(gene_intervals_shuffled_file, 'w') as fileh:
            subprocess.call(args, stdout=fileh)
        # read new gene positions
    #     logging.info("Reading shuffled gene intervals")
        shuffle_gene_map = {}
        with open(gene_intervals_shuffled_file) as fileh:
            for line in fileh:
                print line
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                gene_id = fields[3]
                shuffle_gene_map[gene_id] = (chrom, start, end)
        # reposition transcripts
    #     logging.info("Repositioning transcripts")
        with open(shuffled_bed_file, 'w') as fileh:
            for locus_transcripts in parse_gtf(open(gtf_file)):
                # get original positions
                orig_gene_map = {}
                for g in get_gene_intervals(locus_transcripts):
                    orig_gene_map[g.gene_id] = (g.chrom, g.start, g.end)
                for t in locus_transcripts:
                    gene_id = t.attrs['gene_id']
                    orig_chrom, orig_start, orig_end = orig_gene_map[gene_id]
                    if gene_id not in shuffle_gene_map:
                        logging.warning('Gene %s [%s:%d-%d] could not be shuffled' % (gene_id, orig_chrom, orig_start, orig_end))
                        continue
                    new_chrom, new_start, new_end = shuffle_gene_map[gene_id]
                    # reposition transcript
                    t.chrom = new_chrom
                    t.start = new_start + (t.start - orig_start)
                    t.end = new_start + (t.end - orig_start)
                    for e in t.exons:
                        e.start = new_start + (e.start - orig_start)
                        e.end = new_start + (e.end - orig_start)
                    fields = write_bed(t.chrom, t.attrs['transcript_id'], t.strand, 1000, t.exons)
                    print >>fileh, '\t'.join(map(str,fields))
    #                 for f in t.to_gtf_features(source='shuffle'):
    #                     print >>fileh, str(f)
    
        #gwas snps
        #do intersection for shuffle with GWAS snps
#         logging.info('Performing GWAS intersect %d/%d' % (x+1, NUM_SHUFFS))
        args = ['bedtools', 'intersect', 
                '-a', gwas_bed,
                '-b', shuffled_bed_file,
                '-wa',
                '-wb']
        with open(intersect_file, 'w') as fileh:
            subprocess.call(args, stdout=fileh)
        #count number of GWAS SNPs caught
        snps = set()
        RSIDCOL = 3
        for line in open(intersect_file):
            line = line.strip().split('\t')
            rsID = line[RSIDCOL]
            snps.add(rsID)
        val = len(snps)
        gwas_vals.append(val)
        
        #snp universe
        #do intersections for shuffle and snp universe 
        logging.info('Performing intersects %d/%d' % (x+1, NUM_SHUFFS))
        args = ['bedtools', 'intersect', 
                '-a', snp_bed,
                '-b', shuffled_bed_file,
                '-wa',
                '-wb']
        with open(intersect_file, 'w') as fileh:
            subprocess.call(args, stdout=fileh)
        #count number of SNPs caught
        snps = set()
        RSIDCOL = 3
        for line in open(intersect_file):
            line = line.strip().split('\t')
            rsID = line[RSIDCOL]
            snps.add(rsID)
        snp_overlap = len(snps)
        universe_vals.append(snp_overlap)
        frac = float(val)/snp_overlap
        OR = frac_real/frac
        logging.info('GWAS: %d, Universe: %d, Fraction: %f, OR: %f'%(val, snp_overlap, frac, OR))
        
        print >> valso, '\t'.join(map(str, [val, snp_overlap]))
    valso.close()
    for files in tmp_files: 
        os.remove(files)#         logging.info('Performing GWAS intersect %d/%d' % (x+1, NUM_SHUFFS))


if __name__ == '__main__': 
    sys.exit(main())



