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
import glob
import numpy
from scipy import stats
import plot_test
    
def main():

    
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    
    parser.add_argument("lesion_file")
    parser.add_argument("matrix_file")
    parser.add_argument("output_dir")
    parser.add_argument("tumor_type")
    parser.add_argument("-bed", dest="bed_file", default = '/home/yniknafs/misc_files/compendia_2013_hg18.bed')
    parser.add_argument("-pheno", dest="phenos_file", default = '/home/yniknafs/compendia/phenos.txt')
    parser.add_argument('-RNA', dest="RNA_id_file", default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    parser.add_argument("-narrow", dest='narrow', default = 0)
    parser.add_argument('-at', dest="amp_threshold", default=1.2)
    parser.add_argument('-dt', dest="del_threshold", default=0.9)
    parser.add_argument('-plot', dest="plot", default = 0)
    
    args = parser.parse_args()
    
    folder = os.path.join(args.output_dir, args.tumor_type)
    
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    if not os.path.exists(folder):
        os.mkdir(folder)


#split lesionfile
    amp_lesion_dest = 'amp_lesions.txt'
    del_lesion_dest = 'del_lesions.txt'
    amp_lesion_file = open(amp_lesion_dest, 'w')
    del_lesion_file = open(del_lesion_dest, 'w')
    with open(args.lesion_file, 'r') as f:
        header_fields = f.next().strip().split('\t') 
        reader=csv.reader(f,delimiter='\t')
        amp_lesion_file.write('\t'.join(map(str, header_fields)))
        amp_lesion_file.write('\n')
        del_lesion_file.write('\t'.join(map(str, header_fields)))
        del_lesion_file.write('\n')
        for row in reader: 
            peakID = row[header_fields.index('Unique Name')]
            peak_type = peakID.split(' ')[0]
            if peak_type == 'Amplification':
                amp_lesion_file.write('\t'.join(map(str, row)))
                amp_lesion_file.write('\n')
            if peak_type == "Deletion":
                del_lesion_file.write('\t'.join(map(str,row)))
                del_lesion_file.write('\n')
    amp_lesion_file.close()
    del_lesion_file.close()

    
#generate list of transcripts in each peak
    def bedmaker (lesions, type):
        #make bedfiles to identify overlapping genes from the all_lesions file (in GISTIC output)    
        with open(lesions, 'r') as f:
            labels = f.readline().replace(" ", "").strip().split('\t')
            label_dict = collections.defaultdict(lambda: [])
            for x in xrange(len(labels)): 
                label_dict[labels[x]]=x 
            
            peaks = []
            
            reader=csv.reader(f,delimiter='\t')
            for labels in reader: 
                labels = map(lambda s:s.replace(" ",""), labels)
                if labels[0].endswith("CNvalues") == False:
                    if args.narrow == '1':
                        loc = labels[label_dict['PeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                        peaks.append((loc[0], loc[1], loc[2], labels[label_dict['UniqueName']],labels[label_dict['qvalues']]))
                        print "narrow"
                    else: 
                        loc = labels[label_dict['WidePeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                        peaks.append((loc[0], loc[1], loc[2], labels[label_dict['UniqueName']], labels[label_dict['qvalues']]))
            
            
            file_dest = type + '_tmp.bed'
            
            
            outfile = open(file_dest, 'w')
            for x in peaks: 
                outfile.write('\t'.join(map(str, x)))
                outfile.write('\n')
            outfile.close() 
        return file_dest
    amp_bedfile = bedmaker(amp_lesion_dest, "amp")
    del_bedfile = bedmaker(del_lesion_dest, "del")
    print amp_lesion_dest
    print del_lesion_dest
    
    #perform overlap and create peak file (i.e. list of transcripts in each peak)

    
    def bedtools_intersect(bed_file, type):
        
        peak_file_dest = os.path.join(args.output_dir, type + "_peaks_by_transcript.txt")
        peak_file = open(peak_file_dest, 'w')
        
        p1 = subprocess.check_output(["bedtools", 'intersect', '-a', bed_file, '-b', args.bed_file, '-loj'])
        for line in p1.splitlines():
            items = line.split('\t')
            if items[8] != '.':
                a, b = (items[8].split('|')[1], items[8].split('|')[0])
                outlist = [a, b, items[3], items[0], items[1], items[2]]
                peak_file.write('\t'.join(map(str, outlist)))
                peak_file.write('\n')
        peak_file.close()
        #os.remove(bed_file)
        return peak_file_dest
    
    amp_peak_file_dest = bedtools_intersect(amp_bedfile, 'amp')
    del_peak_file_dest = bedtools_intersect(del_bedfile, 'del')
    
    def id_converter_dict():
    #create DNA_id_file from clinicalMatrix TCGA data (i.e. match TCGA ptID to the CNV sampleID)
    #used to convert CNV_DNA sampleID's to RNA sampleID's which are used in the expression file
        DNA_id_file = open(os.path.join(args.output_dir, "DNA_id_file.txt"), 'w')
        DNA_id_file_dir = os.path.join(args.output_dir, "DNA_id_file.txt")
        with open(args.matrix_file, 'r') as f:
            header_fields = f.next().strip().split('\t') 
            reader=csv.reader(f,delimiter='\t')
            for row in reader: 
                CNV_sampleID = row[header_fields.index('sampleID')]
                ptID = row[header_fields.index('_PATIENT')]
                line = [CNV_sampleID, ptID]
                DNA_id_file.write('\t'.join(map(str, line)))
                DNA_id_file.write('\n')
        DNA_id_file.close()
        
        #building dictionary to convert GISTIC id to library ID   
        pt_to_RNA = collections.defaultdict(lambda: [])
        with open(args.RNA_id_file, 'r') as f:
            next(f)     
            reader=csv.reader(f,delimiter='\t')
            for row in reader: 
                RNA_sample, pt = row
                if pt.startswith("TCGA"):
                    pt_split = pt.split('-')
                    pt_new = '-'.join([pt_split[0], pt_split[1], pt_split[2]])
                    line = [RNA_sample, pt_new]
                    pt_to_RNA[pt_new] = RNA_sample
        DNA_to_RNA = collections.defaultdict(lambda: [])
        with open(DNA_id_file_dir, 'r') as f:
            reader=csv.reader(f,delimiter='\t')
            for row in reader: 
                DNA_sample, pt = row
                DNA_to_RNA[DNA_sample] = pt_to_RNA[pt][:]
        return DNA_to_RNA
    
    DNA_to_RNA_dict = id_converter_dict() 
    
    def id_converter(DNA_id):
        return str(DNA_to_RNA_dict[str(DNA_id)])
    
    
    def amp_id_counter(lesion_file):
               
        test_id_dict = collections.defaultdict(lambda: [])
                
        with open(lesion_file, 'r') as f:
            header_fields = f.next().strip().split('\t') 
            reader=csv.reader(f,delimiter='\t')
            pt_list_all = []
            for x in xrange(9, len(header_fields)):
                if id_converter(header_fields[x]) != '[]':
                    pt_list_all.append(id_converter(header_fields[x]))
            thresh_dict = collections.defaultdict(lambda: [])
            thresh_dict["Amplification"] = (float(args.amp_threshold), 1)
            thresh_dict["Deletion"] = (float(args.del_threshold), -1)
            
            for row in reader:
                
                #print header_fields
                #print header_fields.index('q values')
                #print row[5]
                #print row[header_fields.index('q values')]
                z = header_fields.index('q values')
                #print z
                #print row[z]
                q_value = row[z]
                
                CNV_null = []
                CNV_pos = []
                row.pop()
                stripped_CNV = row[0].strip().replace(" ", "").split('-')
                peakname = stripped_CNV[0]
                stripped_ad = row[0].strip().split(' ')[0]
                if len(stripped_CNV) > 1:
                    for x in xrange(9, len(row)):
                        if id_converter(header_fields[x]) != '[]':
                            product = float(row[x])*float(thresh_dict[stripped_ad][1])
                            if product > thresh_dict[stripped_ad][0]:
                                CNV_pos.append(id_converter(header_fields[x]))
                            else: 
                                CNV_null.append(id_converter(header_fields[x]))
                #q_value = row[header_fields.index('q values')]
                
                test_id_dict[peakname] = [CNV_pos,CNV_null, q_value]
        return test_id_dict, pt_list_all
    amp_id_dict, amp_pt_list_all = amp_id_counter(amp_lesion_dest)
    del_id_dict, del_pt_list_all = amp_id_counter(del_lesion_dest)
# create patient list for generating expression file
    def pt_lister(pt_list, type):
        pt_list_dir = os.path.join(args.output_dir, type + '_pt_list.txt')
        pt_list = open(pt_list_dir, 'w')
        for pt in amp_pt_list_all:
            pt_list.write(pt)
            pt_list.write('\n')
        pt_list.close()
        return pt_list_dir
    amp_pt_dir = pt_lister(amp_pt_list_all, 'amp')
    del_pt_dir = pt_lister(del_pt_list_all, 'del')
        
# create transcript list for generating expression file
    def transcript_lister(peak_file, type):          
        transcript_list_dir = os.path.join(args.output_dir, type + '_transcript_list.txt')
        transcript_list = open(transcript_list_dir, 'w')
        with open(peak_file, 'r') as f:
            reader=csv.reader(f,delimiter='\t')
            for row in reader:
                transcript_list.write(row[0])
                transcript_list.write('\n')
        transcript_list.close()    
        return transcript_list_dir
    amp_transcript_list = transcript_lister(amp_peak_file_dest, 'amp')
    del_transcript_list = transcript_lister(del_peak_file_dest, 'del')
    
# create expression file
    def expression(pt_list, transcript_list, type):
        expression_p = subprocess.Popen(["python", '/home/yniknafs/workspace/CNV/get_isoform_expression.py', pt_list, transcript_list, args.tumor_type + "_" + type])
        expression_p.wait()
        expression_file_dir = os.path.join(args.output_dir, args.tumor_type + "_" + type + '.expr.tsv')
        tsv_file = os.path.join(os.path.dirname(expression_file_dir), os.path.basename(expression_file_dir))
        return tsv_file
    amp_tsv_file = expression(amp_pt_dir, amp_transcript_list, "amp")
    del_tsv_file = expression(del_pt_dir, del_transcript_list, "del")

    
    
#import data for normal tissues/i.e. load phenos file to start
    #write file of normals for tumor type
    normal_ptIDs = open('normal_ptID.txt', 'w')
    normal_ptIDs_dir = os.path.join(args.output_dir, 'normal_ptID.txt')
    
    #counter to count number of normals
    normal_count = 0
    
    with open(args.phenos_file, 'r') as f:
        header_fields = f.next().strip().split('\t')
        reader = csv.reader(f, delimiter = '\t')
        for row in reader: 
            type = row[header_fields.index('tcga_cancer_type')]
            progression = row[header_fields.index('tcga_sample_type')]
            libID = row[header_fields.index('library_id')]

            if type.strip() == args.tumor_type and int(progression) == 11: 
                normal_ptIDs.write(libID)
                normal_ptIDs.write('\n')
                normal_count += 1
    normal_ptIDs.close()    
    amp_norm_tsv_file = expression(normal_ptIDs_dir, amp_transcript_list, 'amp_norm')
    del_norm_tsv_file = expression(normal_ptIDs_dir, del_transcript_list, 'del_norm')
    #generate expression matrix for normals 
 
    #get peak info  
    #make dictionary that converts transcript ID to a peak
    def transcript_to_peak(peak_file):
        transcript_to_peak = collections.defaultdict(lambda: [])
        with open(peak_file, 'r') as f:
            reader=csv.reader(f,delimiter='\t')
            for row in reader: 
                transcriptID = row[0]
                peakID = row[2]
                transcript_to_peak[transcriptID] = peakID
        return transcript_to_peak
    amp_tran_to_peak = transcript_to_peak(amp_peak_file_dest)
    print amp_tran_to_peak
    print amp_peak_file_dest
    del_tran_to_peak = transcript_to_peak(del_peak_file_dest)
    print del_tran_to_peak
    print del_peak_file_dest
    

    
    def tsv_collapser(tsv_file):
        print tsv_file
        gene_expression_file_dir = os.path.join(os.path.dirname(tsv_file), 'gene_'+ os.path.basename(tsv_file))
        print gene_expression_file_dir
        gene_expression_file = open(gene_expression_file_dir, 'w')
        gene_list = collections.defaultdict(lambda: [])
        
        with open(tsv_file, 'r') as f:
            header_fields = f.next().strip().split('\t')
            gene_expression_file.write('\t'.join(map(str, header_fields)))
            gene_expression_file.write('\n')
            
            reader=csv.reader(f,delimiter='\t')        
            i = 0
            gene_row = []
            for row in reader:
                gene = row[header_fields.index('gene_id')]
                if gene not in gene_list.keys():
                    gene_list[gene] = row
                else:
                    working_row = gene_list[gene]
                    for x in xrange(13, len(row)):
                        working_row[x] = float(working_row[x]) + float(row[x])
                        gene_list[gene] = working_row
            for gene in gene_list.iterkeys():
                line = gene_list[gene]
                gene_expression_file.write('\t'.join(map(str, line)))
                gene_expression_file.write('\n')
        gene_expression_file.close()
        return gene_expression_file_dir
    amp_gene_tsv_file = tsv_collapser(amp_tsv_file)
    del_gene_tsv_file = tsv_collapser(del_tsv_file)
    amp_norm_gene_tsv_file = tsv_collapser(amp_norm_tsv_file)
    del_norm_gene_tsv_file = tsv_collapser(del_norm_tsv_file)
    
        
    
    
    def norm_matrix_maker(normal_file_path):
        #create dictionary that provides expression for normals for each transcript
        norm_dict = collections.defaultdict(lambda: [])
        norm_expressions_dict = collections.defaultdict(lambda: [])
        
        with open(normal_file_path, 'r') as f:
            header_fields = f.next().strip().split('\t') 
            reader=csv.reader(f, delimiter='\t')
                
            #make list of all the expression values for the normals
            for row in reader: 
                transcriptID = row[header_fields.index('tracking_id')]
                expressions = []
                for x in xrange(13, len(row)):
                    expressions.append(float(row[x]))
                    
                #mean expression of normal for this transcript
                normal_mean = numpy.mean(expressions)
                #create dictionary entry with the mean expression for each transcript
                norm_expressions_dict[transcriptID] = expressions
        
                #create dictionary that lists all the expression values in order to perform t-test later
                norm_dict[transcriptID] = normal_mean
        return norm_dict, norm_expressions_dict
    
    amp_gene_normal_dict, amp_gene_normal_expressions_dict = norm_matrix_maker(amp_norm_gene_tsv_file)
    del_gene_normal_dict, del_gene_normal_expressions_dict = norm_matrix_maker(del_norm_gene_tsv_file)

    #at this point, for each transcript has a mean expression for normals, and you have the number of normals for each tumor type   
    
    #for each transcript/gene identify average expression in each patient group, perform t-test across two groups, 
    #add in a column for the normal, and perform the fold-change and perform t-test on normal as well

    def tsv_processor(tsv_file, norm_dict, norm_expressions_dict, peak_dict, id_dict, type):
        if args.plot != 0:
            plot_dir = os.path.join(folder, 'plots', type)
            if not os.path.exists(os.path.dirname(plot_dir)):
                os.mkdir(os.path.dirname(plot_dir))
            if not os.path.exists(plot_dir):
                os.mkdir(plot_dir)

        
        file_dir = os.path.join(folder, args.tumor_type + '_' + type + '_processed.cnv')
        cnv_file = open(file_dir, 'w')
        #pos_dict_dir = os.path.join(folder, args.tumor_type + '_' + type + '_dict_pos.txt')
        #null_dict_dir = os.path.join(folder, args.tumor_type + '_' + type + '_dict_null.txt')
        
        with open(tsv_file, 'r') as f:
            header_fields = f.next().strip().split('\t')
            header_fields_copy = header_fields[:13]
            
            if args.plot != 0:
                header_fields_copy[13:13] = ['cancer_type', 'peakID', 'peak_type', 'pos_mean', 'pos_count', 'null_mean', 
                                             'null_count', 'cancer_mean', 'cancer_count', 'normal_mean', 'normal_count',  
                                             'pos_null_fold_change', 'pos_null_p value', 'pos_normal_fold change', 
                                             'pos_normal_p_value', 'cancer_normal_fold_change', 'cancer_normal_p_value',
                                             'plot_hyperlink', 'plot_log_hyperlink', 'hyperlink_antecedent', 'log_hyperlink_antecedent'] 
            else:
                header_fields_copy[13:13] = ['cancer_type', 'peakID', 'peak_type', 'pos_mean', 'pos_count', 'null_mean', 
                                             'null_count', 'cancer_mean', 'cancer_count', 'normal_mean', 'normal_count',  
                                             'pos_null_fold_change', 'pos_null_p value', 'pos_normal_fold change', 
                                             'pos_normal_p_value', 'cancer_normal_fold_change', 'cancer_normal_p_value'] 
            cnv_file.write('\t'.join(map(str, header_fields_copy)))
            cnv_file.write('\n')
            
            
            
            reader=csv.reader(f,delimiter='\t')        
            i = 2
            protein_count = collections.defaultdict(lambda:0)
            ncrna_count = collections.defaultdict(lambda:0)
            genes_in_peak = collections.defaultdict(lambda: [])
            
            for row in reader:
                print "tit"
                pos_expr = []
                null_expr = []
                transcript = row[0]
                peak = peak_dict[transcript]   
                
                
                #generating count for number of pts with CNV or number of pts w/o CNV
                #don't forget we have the variable "norm_counter" above
                print type
                print id_dict
                print peak
                
                pos_count = len(id_dict[peak][0])
                null_count = len(id_dict[peak][1])
                cancer_count = pos_count + null_count
                
                ssea_ids = []
                ssea_weights = []
                
                #making list of expression values for each pt with the CNV for this transcript
                for x in id_dict[peak][0]:
                    pos_expr.append(float(row[header_fields.index(x)]))
                    ssea_ids.append(x)
                    ssea_weights.append(float(row[header_fields.index(x)]))
                #making list of expression values for each pt without the CNV for this transcript
                for x in id_dict[peak][1]:
                    null_expr.append(float(row[header_fields.index(x)]))
                    ssea_ids.append(x)
                    ssea_weights.append(float(row[header_fields.index(x)]))
                
                              
                #means 
                pos_mean = numpy.mean(pos_expr)
                null_mean = numpy.mean(null_expr)
                normal_mean = norm_dict[transcript]
                print norm_dict
                cancer_mean = (pos_mean*pos_count + null_mean*null_count)/cancer_count
                
                
                #calculate fold changes
                print pos_mean
                print null_mean
                pos_null_fc = (pos_mean+0.0000000001)/(null_mean+0.0000000001)
                print pos_mean
                #print normal_mean
                print transcript
                pos_normal_fc = (pos_mean+0.0000000001)/(normal_mean+0.0000000001)
                cancer_normal_fc = (cancer_mean+0.0000000001)/(normal_mean+0.0000000001)
                
                
                #generate lists for t-test 
                cancer_expr = pos_expr + null_expr
                normal_expr = norm_expressions_dict[transcript]
                
                #perform t-tests
                t_pos_null,p_pos_null = stats.ttest_ind(pos_expr, null_expr)
                t_pos_normal,p_pos_normal = stats.ttest_ind(pos_expr, normal_expr)
                t_cancer_normal,p_cancer_normal = stats.ttest_ind(cancer_expr, normal_expr)      
                #p_pos_null = 4
                #p_pos_normal = 4
                #p_cancer_normal = 4
                #create peak-type
                peak_type = peak.split('Peak')[0]
                
                #lists of expression values for each category
                #normal --> norm_expressions_dict[transcript]
                #pos --> pos_expr
                #null --> null_expr
                
                if args.plot != 0:
                    plot_file_dir = os.path.abspath(os.path.join(plot_dir, transcript + '.png'))
                    plot_file_log_dir = os.path.abspath(os.path.join(plot_dir, transcript + '_log.png'))
                    
                    plot_test.plotter(pos_expr, null_expr, norm_expressions_dict[transcript], transcript, plot_file_dir)
                    plot_hyperlink = '=HYPERLINK("file://"&AG' + str(i) + '; \"non-log\"'
                    plot_test.plotter(pos_expr, null_expr, norm_expressions_dict[transcript], transcript + " log", plot_file_log_dir, True)
                    plot_log_hyperlink = '=HYPERLINK("file://"&AH' + str(i) + '; \"log\"'
                    
                    
                    row_clip = row[:13]
                    row_clip[13:13] = [args.tumor_type, peak, peak_type, pos_mean, pos_count, null_mean, null_count, cancer_mean, cancer_count, 
                                  normal_mean, normal_count, pos_null_fc, p_pos_null, pos_normal_fc, p_pos_normal,
                                  cancer_normal_fc, p_cancer_normal, plot_hyperlink, plot_log_hyperlink, plot_file_dir, plot_file_log_dir ]
                    cnv_file.write('\t'.join(map(str, row_clip)))        
                    cnv_file.write("\n")
                else:
                    row_clip = row[:13]
                    row_clip[13:13] = [args.tumor_type, peak, peak_type, pos_mean, pos_count, null_mean, null_count, cancer_mean, cancer_count, 
                                  normal_mean, normal_count, pos_null_fc, p_pos_null, pos_normal_fc, p_pos_normal,
                                  cancer_normal_fc, p_cancer_normal]
                    cnv_file.write('\t'.join(map(str, row_clip)))        
                    cnv_file.write("\n")
                
                
                
                category = row[header_fields.index('category')]
                if category == "protein":
                    protein_count[peak]+=1
                else:
                    ncrna_count[peak]+=1
                gene = row[header_fields.index('gene_id')]
                genes_in_peak[peak].append(gene)
                i+=1
        cnv_file.close()
        return ncrna_count, protein_count, genes_in_peak

    amp_ncrna_count_dict, amp_protein_count_dict, amp_genes_in_peak_dict = tsv_processor(amp_gene_tsv_file, 
        amp_gene_normal_dict, amp_gene_normal_expressions_dict, 
        amp_tran_to_peak, amp_id_dict, "amp")

    del_ncrna_count_dict, del_protein_count_dict, del_genes_in_peak_dict = tsv_processor(del_gene_tsv_file, 
        del_gene_normal_dict, del_gene_normal_expressions_dict, 
        del_tran_to_peak, del_id_dict, "del")

    def peak_processor(lesion_file, ncrna_count_dict, protein_count_dict, genes_in_peak_dict, type):
        peak_cnv_file_dir = os.path.join(folder, args.tumor_type + '_' + type +'_peaks_processed.cnv')
        peak_cnv_file = open(peak_cnv_file_dir, 'w')                
        with open(lesion_file, 'r') as f:
            header_fields = f.next().strip().split('\t') 
            header_fields_copy = header_fields[:7]
            header_fields_copy[7:7] = ['proteins_in_peak', 'ncrnas_in_peak', 'genes_in_peak']
            peak_cnv_file.write('\t'.join(map(str,header_fields_copy)))
            peak_cnv_file.write('\n') 
            reader=csv.reader(f,delimiter='\t')
            
            for row in reader:
                CNV_null = []
                CNV_pos = []
                row.pop()
                stripped_CNV = row[0].strip().replace(" ", "").split('-')
                peakname = stripped_CNV[0]
                stripped_ad = row[0].strip().split(' ')[0]
                ncrna_count = ncrna_count_dict[peakname]
                protein_count = protein_count_dict[peakname]
                genes_in_peak = genes_in_peak_dict[peakname]
                if len(stripped_CNV) > 1:
                    row_clip = row[:7]
                    row_clip[7:7] = [protein_count, ncrna_count, ','.join(genes_in_peak)]
                    
                    peak_cnv_file.write('\t'.join(map(str,row_clip)))
                    peak_cnv_file.write('\n')
        peak_cnv_file.close()
        
    peak_processor(amp_lesion_dest, amp_ncrna_count_dict, amp_protein_count_dict, amp_genes_in_peak_dict, "amp")
    peak_processor(del_lesion_dest, del_ncrna_count_dict, del_protein_count_dict, del_genes_in_peak_dict, "del")

#     
    #p = subprocess.call(['rm', peak_file_dest, peak_file_other_dest, DNA_id_file_dir, lesion_processed_dir, pt_list_dir, transcript_list_dir, normal_ptIDs_dir])    
    #p = subprocess.call('rm *.tsv', shell=True)
#     

    
    
    
    
    logging.info("Finished")   
    


    return 0



if __name__ == '__main__':

    sys.exit(main())