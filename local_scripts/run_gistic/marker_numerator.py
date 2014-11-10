import os
import sys
import csv
import time
import collections
import bisect






fname_probes = "D:\yniknafs\CNV\\TCGA.DCC.GenomeWideSNP6.marker.na30.txt"
file_direct = "D:\yniknafs\CNV\Real_samples\CNV_raw"
filelist = os.listdir(file_direct)


def numerator(probes_file, CNV_file):
    probes = collections.defaultdict(lambda: [])
    
    with open(probes_file) as f:
        next(f)
        reader=csv.reader(f,delimiter='\t')
        for name, chrs, position in reader:
            probes[chrs].append(int(position))
  
            
    for chrom in probes.iterkeys():
        probes[chrom] = sorted(probes[chrom])
    
    newsegfile = CNV_file[:-4] + "_withnummarkers.txt"
    outfile = open(newsegfile, 'w')
    
    with open(CNV_file) as f:
        reader=csv.reader(f,delimiter='\t')
        for name, chrs, start, end, strand, score in reader:
            chrs = chrs[3:]
            start = int(start)
            end = int(end)
            cnv_probes = probes[chrs]
            start_index = bisect.bisect_left(cnv_probes, start)
            end_index = bisect.bisect_left(cnv_probes, end)
            numprobes = end_index - start_index
            
            line = name, chrs, start, end, numprobes, score
            outfile.write('\t'.join(map(str, line)))
            outfile.write('\n')
            
    outfile.close()


for file in filelist: 
    filename = file_direct + "\\" + file 
    numerator(fname_probes, filename)