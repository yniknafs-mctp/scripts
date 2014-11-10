import argparse
import logging
import os
import sys
import csv
import time
import collections
import bisect
import operator
import subprocess
import numpy
from scipy import stats

def main():

    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    logging.info("Running code")
    parser = argparse.ArgumentParser()
    parser.add_argument("transcript_file")
    parser.add_argument("peak_file")
    parser.add_argument("expression_file")
    args = parser.parse_args()
    
    #building dictionary to convert GISTIC id to library ID   
    transcript_peak = collections.defaultdict(lambda: [])
    peak_info = collections.defaultdict(lambda: [])
    with open(args.transcript_file, 'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            transcriptID = row[0]
            peakID = row[2]
            transcript_peak[transcriptID] = peakID
                
    with open(args.peak_file, 'r') as f:
        header_fields = f.next().strip().split('\t') 
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            peakID = row[header_fields.index('Unique Name')]
            peakID_trimmed = peakID.strip().replace(' ',''). replace ('-CNvalues','')
            CNV_null_c_handle = row[header_fields.index('CNV_null_c')]
            CNV_pos_c_handle = row[header_fields.index('CNV_pos_c')]
            CNV_null_c = CNV_null_c_handle.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
            CNV_pos_c = CNV_pos_c_handle.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
            peak_lim = row[header_fields.index('Peak Limits')]
            q_val = row[header_fields.index('q values')]
            peak_info[peakID_trimmed] = [CNV_pos_c, CNV_null_c, peak_lim, q_val]
    
    
    outfile1 = open('transcript_processed.txt', 'w')
    outfile2 = open('genes_processed.txt', 'w')
    with open(args.expression_file, 'r') as f:
        header_fields = f.next().strip().split('\t')
        header_fields_copy = header_fields[:]
        header_fields_copy[13:13] = ['peakID', 'pos_mean', 'null_mean', 'fold_change', 'p value'] 
        outfile1.write('\t'.join(map(str, header_fields_copy)))
        outfile1.write('\n')
        outfile2.write('\t'.join(map(str, header_fields_copy)))
        outfile2.write('\n')
        reader=csv.reader(f,delimiter='\t')
        i = 0
        gene_row = []
        for row in reader:
            if i == 0: 
                gene_row = row[:]
            if i != 0:
                if row[6] == gene_row[6]:
                    for x in xrange(13, len(row)):
                        gene_row[x] = float(gene_row[x]) + float(row[x])
                        print row[x]
                else:
                    pos_expr_gene = []
                    null_expr_gene = []
                    transcript_gene = gene_row[0]
                    peak_gene = transcript_peak[transcript_gene]
                    for x in peak_info[peak_gene][0]:
                        pos_expr_gene.append(float(gene_row[header_fields.index(x)]))
                    for x in peak_info[peak_gene][1]:
                        null_expr_gene.append(float(gene_row[header_fields.index(x)]))
                    pos_mean_gene = numpy.mean(pos_expr_gene)
                    null_mean_gene = numpy.mean(null_expr_gene)
                    fold_change_gene = (pos_mean_gene+0.0000000001)/(null_mean_gene+0.0000000001)
                    t_gene,p_gene = stats.ttest_ind(pos_expr_gene, null_expr_gene)
                    gene_row[13:13] = [peak_gene, pos_mean_gene, null_mean_gene, fold_change_gene, p_gene]
                    outfile2.write('\t'.join(map(str, gene_row)))
                    outfile2.write('\n')
                    gene_row = row[:] 
                    
                
            

            pos_expr = []
            null_expr = []
            transcript = row[0]
            peak = transcript_peak[transcript]
            for x in peak_info[peak][0]:
                pos_expr.append(float(row[header_fields.index(x)]))
            for x in peak_info[peak][1]:
                null_expr.append(float(row[header_fields.index(x)]))
            pos_mean = numpy.mean(pos_expr)
            null_mean = numpy.mean(null_expr)
            fold_change = (pos_mean+0.0000000001)/(null_mean+0.0000000001)
            t,p = stats.ttest_ind(pos_expr, null_expr)
            row[13:13] = [peak, pos_mean, null_mean, fold_change, p]
            outfile1.write('\t'.join(map(str, row)))        
            outfile1.write("\n")
            i+=1
        pos_expr_gene = []
        null_expr_gene = []
        transcript_gene = gene_row[0]
        peak_gene = transcript_peak[transcript_gene]
        for x in peak_info[peak_gene][0]:
            pos_expr_gene.append(float(gene_row[header_fields.index(x)]))
        for x in peak_info[peak_gene][1]:
            null_expr_gene.append(float(gene_row[header_fields.index(x)]))
        pos_mean_gene = numpy.mean(pos_expr_gene)
        null_mean_gene = numpy.mean(null_expr_gene)
        fold_change_gene = (pos_mean_gene+0.0000000001)/(null_mean_gene+0.0000000001)
        t_gene,p_gene = stats.ttest_ind(pos_expr_gene, null_expr_gene)
        gene_row[13:13] = [peak_gene, pos_mean_gene, null_mean_gene, fold_change_gene, p_gene]
        outfile2.write('\t'.join(map(str, gene_row)))
        outfile2.write('\n')
    outfile1.close()
    outfile2.close()
#         header_fields = f.next().strip().split('\t') 
#         reader=csv.reader(f,delimiter='\t')
#         
#         header_fields_copy = header_fields[:]
#         header_fields_copy[9:9] = ["CNV_pos", "CNV_pos_c", "CNV_null_c", "CNV_null"]
#         outfile.write('\t'.join(map(str, header_fields_copy)))
#         outfile.write('\n')
#         thresh_dict = collections.defaultdict(lambda: [])
#         thresh_dict["Amplification"] = (args.amp_threshold, 1)
#         thresh_dict["Deletion"] = (args.del_threshold, -1)
#         for row in reader:
#             CNV_null = []
#             CNV_pos = []
#             CNV_null_c = []
#             CNV_pos_c = []
#             row.pop()
#             stripped_CNV = row[0].strip().replace(" ", "").split('-')
#             stripped_ad = row[0].strip().split(' ')[0]
#             if len(stripped_CNV) > 1:
#                 for x in xrange(9, len(row)):
#                     product = float(row[x])*float(thresh_dict[stripped_ad][1])
#                     if product > thresh_dict[stripped_ad][0]:
#                         CNV_pos.append(header_fields[x])
#                         if DNA_to_RNA[str(header_fields[x])] != []:
#                             CNV_pos_c.append(DNA_to_RNA[str(header_fields[x])])
#                     else: 
#                         CNV_null.append(header_fields[x])
#                         if DNA_to_RNA[str(header_fields[x])] != []:
#                             CNV_null_c.append(DNA_to_RNA[str(header_fields[x])])
#                 row[9:9] = [CNV_pos, CNV_pos_c, CNV_null, CNV_null_c]
#                 outfile.write('\t'.join(map(str, row)))        
#                 outfile.write("\n")
#         
#     outfile.close()
# #             
    
    
    logging.info("Finished")   

    return 0



if __name__ == '__main__':

    sys.exit(main())