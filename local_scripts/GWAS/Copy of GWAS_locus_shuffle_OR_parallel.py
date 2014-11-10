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
1) will shuffle the intergenic LOCI inside the intergenic space (compared to ref merged bed)
2) will intersect shuffle with gwas snps to check for overlap
3) for each shuffle report number of GWAS snps that are caught, and OR for real
compared to shuffle for frac of GWAS snps caught/All snps caught

does this in parallel with multi processor support
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


def shuffle(process,
            snps_file,
            gwas_file,
            excl_file2, 
            locus_intervals_file, 
            chrom_sizes_file,
            gtf_file2,
            frac_real,
            gwas_real,
            snps_real,
            NUM_SHUFFS,
            output_dir,
            flank,
            exon):
    x = process
    prefix = 'process' + str(x)
    locus_intervals_shuffled_file = os.path.join(output_dir, prefix + '.locus_intervals.shuffle.bed')
    shuffled_bed_file = os.path.join(output_dir, prefix + '.shuffle.bed')
    intersect_file = os.path.join(output_dir, prefix + '.intersect.txt')
    args_shuff = ['bedtools', 'shuffle', 
            '-excl', excl_file2,
            '-i', locus_intervals_file, 
            '-g', chrom_sizes_file]
    with open(locus_intervals_shuffled_file, 'w') as fileh:
        subprocess.call(args_shuff, stdout=fileh)
    # read new gene positions
#     logging.info("Reading shuffled gene intervals")
    shuffle_locus_map = {}
    with open(locus_intervals_shuffled_file) as fileh:
        for line in fileh:
            fields = line.strip().split('\t')
#             print fields
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            locus_id = int(fields[3])
            shuffle_locus_map[locus_id] = (chrom, start, end)
    # reposition transcripts
    
    #read chrom file to make sure flanks added do not enter chrom ends
    chrom_length = {}
    for line in open(chrom_sizes_file): 
        line = line.strip().split('\t')
        chr = line[0]
        length = line[1]
        chrom_length[chr] = length
    
    with open(shuffled_bed_file, 'w') as fileh:
        i=0
        for locus_transcripts in parse_gtf(open(gtf_file2)):
            locus_chrom = locus_transcripts[0].chrom
            locus_start = min(t.start for t in locus_transcripts)
            locus_end = max(t.end for t in locus_transcripts)
            locus_id = i
            i+=1
            orig_locus_map = {}
            orig_locus_map[locus_id] = (locus_chrom, locus_start, locus_end)
            orig_chrom, orig_start, orig_end = orig_locus_map[locus_id]
            if locus_id not in shuffle_locus_map.keys():
                    logging.warning('Locus %s [%s:%d-%d] could not be shuffled' % (locus_id, orig_chrom, orig_start, orig_end))
                    continue
            
            for t in locus_transcripts:
                new_chrom, new_start, new_end = shuffle_locus_map[locus_id]
                # reposition transcript
                t.chrom = new_chrom
                t.start = new_start + (t.start - orig_start)
                t.end = new_start + (t.end - orig_start)
                for e in t.exons:
                    e.start = new_start + (e.start - orig_start)
                    e.end = new_start + (e.end - orig_start)
                fields = write_bed(t.chrom, t.attrs['transcript_id'], t.strand, 1000, t.exons, flank, chrom_length)
                print >>fileh, '\t'.join(map(str,fields))


    #gwas snps
    #do intersection for shuffle with GWAS snps
#         logging.info('Performing GWAS intersect %d/%d' % (x+1, NUM_SHUFFS))
    args_int = ['bedtools', 'intersect', 
            '-a', gwas_file,
            '-b', shuffled_bed_file,
            '-wa',
            '-wb']
    if exon:
        args_int.append('-split')
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args_int, stdout=fileh)
    #count number of GWAS SNPs caught
    snps = set()
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    val = len(snps)
    
    #snp universe
    #do intersections for shuffle and snp universe 
    args_int = ['bedtools', 'intersect', 
            '-a', snps_file,
            '-b', shuffled_bed_file,
            '-wa',
            '-wb']
    if exon:
        args_int.append('-split')
    with open(intersect_file, 'w') as fileh:
        subprocess.call(args_int, stdout=fileh)
    #count number of SNPs caught
    snps = set()
    for line in open(intersect_file):
        line = line.strip().split('\t')
        rsID = line[RSIDCOL]
        snps.add(rsID)
    snp_overlap = len(snps)
    frac = float(val)/snp_overlap
    OR = frac_real/frac
    logging.info('Shuffle %d/%d. GWAS: %d, Universe: %d, Fraction: %f, OR: %f'%(x, NUM_SHUFFS, val, snp_overlap, frac, OR))
    return '\t'.join(map(str, [val, snp_overlap, frac, gwas_real, snps_real, frac_real, OR]))

def shuffle_imap(argtuple):
    return shuffle(*argtuple)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--rand_snp", dest = 'rand_snp',
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
    parser.add_argument("--shuffs", dest = 'shuffs',
                    default = 100,
                    help = 'number of shuffles to perform')
    parser.add_argument("-p", dest = 'proc',
                    default = 4,
                    help = 'number of processors to use')
    parser.add_argument("--flank", dest = 'flank',
                    default = 0,
                    help = 'number of flanking bases to add to bed files')
    parser.add_argument("--exon", dest="exon", 
                        action="store_true", default=False, 
                        help="Perform analysis looking only at exonic overlap")
    args = parser.parse_args()
    
    args.proc = int(args.proc)
    args.flank = int(args.flank)
    
    logging.info('Output is printed to stdout')
    if args.exon: 
        logging.info('Looking at exonic overlap only')
    # check command line parameters
    if which('bedtools') is None:
        parser.error('bedtools binary not found in PATH')
    if not os.path.exists(chrom_sizes_file):
        parser.error('chrom sizes file %s not found' % (chrom_sizes_file))

    if not os.path.isdir('GWAS_TMPS'):
        os.mkdir('GWAS_TMPS')
    
    prefix = 'GWAS_TMPS'    
    locus_intervals_file = os.path.join(prefix, 'locus_intervals.bed')
    intersect_file = os.path.join(prefix, 'intersect.txt')
    assembly_flank = os.path.join(prefix, 'flank.bed')
    assembly_bed = os.path.join(prefix, 'assembly.bed')
    
    #read chrom file to make sure flanks added do not enter chrom ends
    chrom_length = {}
    for line in open(chrom_sizes_file): 
        line = line.strip().split('\t')
        chr = line[0]
        length = line[1]
        chrom_length[chr] = length
    
    
    #convert GTF file to BED for initial intersections 
    logging.info('Converting gtf to BED for initial intersection')
    with open(assembly_bed, 'w') as f:
        for transcripts in parse_gtf(open(args.gtf)):
            for t in transcripts:
                name = t.attrs['transcript_id']
                fields = write_bed(t.chrom, name, t.strand, 1000, t.exons, args.flank, chrom_length)
                print >>f, '\t'.join(fields)
        
    logging.info('Parsing GTF file')
    with open(locus_intervals_file, 'w') as f:
        j = 0
        for locus_transcripts in parse_gtf(open(args.gtf)):
            # find borders of locus
            locus_chrom = locus_transcripts[0].chrom
            locus_start = min(t.start for t in locus_transcripts)
            locus_end = max(t.end for t in locus_transcripts)
            locus_id = j
            j+=1
            print >>f, '\t'.join(map(str, [locus_chrom, locus_start, locus_end, locus_id]))   
    #apply flank to the bed file 
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
    if args.exon:
        args_int.append('-split')
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
            '-wb']
    if args.exon:
        args_int.append('-split')
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
    
    
    #loop the shuffle to generate a distribution of nulls for number of snps hit by random intergenic genes
    pool = multiprocessing.Pool(args.proc)
    NUM_SHUFFS = int(args.shuffs)
    shuff_args = (args.snps,
                  args.gwas,
                  args.excl,
                  locus_intervals_file,
                  args.chrom,
                  args.gtf,
                  frac_real,
                  gwas_overlap,
                  snp_overlap,
                  NUM_SHUFFS,
                  prefix,
                  args.flank,
                  args.exon)
    tasks = []
    header = [
              'gwas_shuff',
              'snp_shuff',
              'frac_shuff',
              'gwas_real',
              'snp_real',
              'frac_real',
              'OR'
              ]
    print '\t'.join(header)
    for i in xrange(NUM_SHUFFS):
        tasks.append((i,) + shuff_args)
    result_iter = pool.imap_unordered(shuffle_imap, tasks)
    for line in result_iter:
        print line
    pool.close()
    pool.join()
    
#     shutil.rmtree(prefix)
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())



