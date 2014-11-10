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
    parser.add_argument("-tl", dest="transcript_list")
    parser.add_argument("-toi", dest="toi")
    parser.add_argument("-bed", dest="bed_file", default = '/home/yniknafs/misc_files/compendia_2013_hg18.bed')
    parser.add_argument("-pheno", dest="phenos_file", default = '/home/yniknafs/compendia/phenos.txt')
    parser.add_argument('-RNA', dest="RNA_id_file", default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    parser.add_argument("-narrow", dest='narrow', default = 1)
    parser.add_argument("output_dir")
    parser.add_argument('-at', dest="amp_threshold", default=0.1)
    parser.add_argument('-dt', dest="del_threshold", default=0.1)
    parser.add_argument('-name', dest="tumor_type")
    args = parser.parse_args()
    
    folder = os.path.join(args.output_dir, args.tumor_type)
    
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    if not os.path.exists(folder):
        os.mkdir(folder)

    
#generate list of transcripts in each peak
    
    #make bedfiles to identify overlapping genes    
    with open(args.lesion_file, 'r') as f:
        labels = f.readline().replace(" ", "").strip().split('\t')
        label_dict = collections.defaultdict(lambda: [])
        for x in xrange(len(labels)): 
            label_dict[labels[x]]=x 
        
        narrow_peaks = []
        wide_peaks = []
        
        reader=csv.reader(f,delimiter='\t')
        for labels in reader: 
            labels = map(lambda s:s.replace(" ",""), labels)
            if labels[0].endswith("CNvalues") == False:
                np_loc = labels[label_dict['PeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                wp_loc = labels[label_dict['WidePeakLimits']].replace('(',':').replace(')', '').replace('-',':').split(':')
                narrow_peaks.append((np_loc[0], np_loc[1], np_loc[2], labels[label_dict['UniqueName']],labels[label_dict['qvalues']]))
                wide_peaks.append((wp_loc[0], wp_loc[1], wp_loc[2], labels[label_dict['UniqueName']], labels[label_dict['qvalues']]))
        outfile_n = open('narrow.tmp.bed', 'w')
        outfile_w = open('wide.tmp.bed', 'w')
        for x in narrow_peaks: 
            outfile_n.write('\t'.join(map(str, x)))
            outfile_n.write('\n')
        for x in wide_peaks: 
            outfile_w.write('\t'.join(map(str, x)))
            outfile_w.write('\n')
        outfile_n.close()
        outfile_w.close() 

    
    

    
    #perform overlap and create peak file (i.e. list of transcripts in each peak)
    peak_file_narrow = open(os.path.join(args.output_dir, "narrowpeaks_by_transcript.txt"), 'w')
    peak_file_wide = open(os.path.join(args.output_dir, "widepeaks_by_transcript.txt"), 'w')
    
    p1 = subprocess.check_output(["bedtools", 'intersect', '-a', 'narrow.tmp.bed', '-b', args.bed_file, '-loj'])
    p2 = subprocess.check_output(["bedtools", 'intersect', '-a', 'wide.tmp.bed', '-b', args.bed_file, '-loj'])
    
    for line in p1.splitlines():
        items = line.split('\t')
        if items[8] != '.':
            a, b = (items[8].split('|')[1], items[8].split('|')[0])
            outlist = [a, b, items[3], items[0], items[1], items[2]]
            peak_file_narrow.write('\t'.join(map(str, outlist)))
            peak_file_narrow.write('\n')
    peak_file_narrow.close()
    
    for line in p2.splitlines():
        items = line.split('\t')
        if items[8] != '.':
            a, b = (items[8].split('|')[1], items[8].split('|')[0])
            outlist = [a, b, items[3], items[0], items[1], items[2]]
            peak_file_wide.write('\t'.join(map(str, outlist)))
            peak_file_wide.write('\n')
    peak_file_wide.close()
    
    os.remove("narrow.tmp.bed")
    os.remove("wide.tmp.bed")
    
    
    #identify which transcript file to use for further processing (i.e. use narrow or wide peaks)
    if args.narrow == 1:
        peak_file_dest = os.path.join(args.output_dir, "narrowpeaks_by_transcript.txt")
        peak_file_other_dest = os.path.join(args.output_dir, "widepeaks_by_transcript.txt") 
    else:
        peak_file_dest = os.path.join(args.output_dir, "widepeaks_by_transcript.txt")
        peak_file_other_dest = os.path.join(args.output_dir, "narrowpeaks_by_transcript.txt")

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


#process all_lesion file to provide list of patients with/without CNV at each peak
    
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
            
    
    #creating processed file
    lesion_processed = open(os.path.join(args.output_dir, 'lesion_processed.txt'), 'w')
    lesion_processed_dir = os.path.join(args.output_dir, 'lesion_processed.txt')
    
    with open(args.lesion_file, 'r') as f:
        header_fields = f.next().strip().split('\t') 
        reader=csv.reader(f,delimiter='\t')
        
        header_fields_copy = header_fields[:]
        header_fields_copy[9:9] = ["CNV_pos", "CNV_pos_c","CNV_null","CNV_null_c"]
        lesion_processed.write('\t'.join(map(str, header_fields_copy)))
        lesion_processed.write('\n')
        thresh_dict = collections.defaultdict(lambda: [])
        thresh_dict["Amplification"] = (float(args.amp_threshold), 1)
        thresh_dict["Deletion"] = (float(args.del_threshold), -1)
        for row in reader:
            CNV_null = []
            CNV_pos = []
            CNV_null_c = []
            CNV_pos_c = []
            row.pop()
            stripped_CNV = row[0].strip().replace(" ", "").split('-')
            stripped_ad = row[0].strip().split(' ')[0]
            if len(stripped_CNV) > 1:
                for x in xrange(9, len(row)):
                    product = float(row[x])*float(thresh_dict[stripped_ad][1])
                    if product > thresh_dict[stripped_ad][0]:
                        CNV_pos.append(header_fields[x])
                        if DNA_to_RNA[str(header_fields[x])] != []:
                            CNV_pos_c.append(DNA_to_RNA[str(header_fields[x])])
                    else: 
                        CNV_null.append(header_fields[x])
                        if DNA_to_RNA[str(header_fields[x])] != []:
                            CNV_null_c.append(DNA_to_RNA[str(header_fields[x])])
                row[9:9] = [CNV_pos, CNV_pos_c, CNV_null, CNV_null_c]
                lesion_processed.write('\t'.join(map(str, row)))        
                lesion_processed.write("\n")
    lesion_processed.close()
    

# create patient list for generating expression file
    pt_list = open(os.path.join(args.output_dir, 'pt_list.txt'), 'w')
    pt_list_dir = os.path.join(args.output_dir, 'pt_list.txt')
    with open(lesion_processed_dir, 'r') as f:
        next(f)     
        reader=csv.reader(f,delimiter='\t')
        a = next(f).strip().split("\t") 
        newlist = a[10][:]
        newlist2 = newlist.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
        for x in newlist2: 
            pt_list.write(x)
            pt_list.write('\n')
        newlist3 = a[12][:]
        newlist4 = newlist3.replace('[', '').replace(']','').replace(',','\t').replace('\'','').strip().replace(' ','').split('\t')
        for x in newlist4:
            pt_list.write(x)
            pt_list.write('\n')
    pt_list.close()
    
          

# create expression file
    expression_p = subprocess.Popen(["python", '/home/yniknafs/workspace/CNV/get_isoform_expression.py', pt_list_dir, args.transcript_list, args.tumor_type])
    expression_p.wait()


# generate final transcript_processed/gene_processed files
    
    #get peak info  
    transcript_peak = collections.defaultdict(lambda: [])
    peak_info = collections.defaultdict(lambda: [])
    with open(peak_file_dest, 'r') as f:
        reader=csv.reader(f,delimiter='\t')
        for row in reader: 
            transcriptID = row[0]
            peakID = row[2]
            transcript_peak[transcriptID] = peakID

    #get list of patients (CNV or normal) for each peak and list some peak info          
    with open(lesion_processed_dir, 'r') as f:
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
    
    #generate expression matrix for normals 
    
    in_var =  'normal_' + args.tumor_type
    expression_p2 = subprocess.Popen(["python", '/home/yniknafs/workspace/CNV/get_isoform_expression.py', normal_ptIDs_dir, args.transcript_list, in_var])
    expression_p2.wait()
    
    
    def tsv_collapser(tsv_file):
        gene_expression_file_dir = os.path.join(os.path.dirname(tsv_file), 'gene_'+ os.path.basename(tsv_file))
        gene_expression_file = open(gene_expression_file_dir, 'w')
        
        
        with open(tsv_file, 'r') as f:
            header_fields = f.next().strip().split('\t')
            gene_expression_file.write('\t'.join(map(str, header_fields)))
            gene_expression_file.write('\n')
            
            reader=csv.reader(f,delimiter='\t')        
            i = 0
            gene_row = []
            for row in reader:
                if i == 0: 
                    gene_row = row[:]
                else:
                    if row[6] == gene_row[6]:
                        for x in xrange(13, len(row)):
                            gene_row[x] = float(gene_row[x]) + float(row[x])
                    else:
                        gene_expression_file.write('\t'.join(map(str, gene_row)))
                        gene_expression_file.write('\n')
                        gene_row = row[:] 
                i+=1         
            gene_expression_file.write('\t'.join(map(str, gene_row)))
            gene_expression_file.write('\n')
        gene_expression_file.close()
        
    normal_expression_file_dir = os.path.join(args.output_dir, 'normal_' + args.tumor_type + '.expr.tsv')
    gene_normal_expression_file_dir = os.path.join(os.path.dirname(normal_expression_file_dir), 'gene_'+ os.path.basename(normal_expression_file_dir))
    tsv_collapser(normal_expression_file_dir)
    
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
    
    normal_dict, normal_expressions_dict = norm_matrix_maker(normal_expression_file_dir)
    gene_normal_dict, gene_normal_expressions_dict = norm_matrix_maker(gene_normal_expression_file_dir)

    #at this point, for each transcript has a mean expression for normals, and you have the number of normals for each tumor type   


    #create a gene_file that is the collapsed .tsv transcript file 
    expression_file_dir = os.path.join(args.output_dir, args.tumor_type + '.expr.tsv')
    
    tsv_collapser(expression_file_dir)
    
    
    #for each transcript/gene identify average expression in each patient group, perform t-test across two groups, 
    #add in a column for the normal, and perform the fold-change and perform t-test on normal as well

    def tsv_processor(type):
        plot_dir = os.path.join(folder, 'nearby_plots', type)
        
        if not os.path.exists(os.path.dirname(plot_dir)):
            os.mkdir(os.path.dirname(plot_dir))
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)

        if type == 'gene':
            tsv_file = os.path.join(os.path.dirname(expression_file_dir), 'gene_'+ os.path.basename(expression_file_dir))
            norm_dict = gene_normal_dict
            norm_expressions_dict = gene_normal_expressions_dict
        elif type == 'transcript': 
            tsv_file = expression_file_dir
            norm_dict = normal_dict
            norm_expressions_dict = normal_expressions_dict
        else: 
            print "poo"
        
        file_dir = os.path.join(folder, args.tumor_type + '_' + type + '_processed.cnv')
        #pos_dict_dir = os.path.join(folder, args.tumor_type + '_' + type + '_dict_pos.txt')
        #null_dict_dir = os.path.join(folder, args.tumor_type + '_' + type + '_dict_null.txt')
        
        cnv_file = open(file_dir, 'w')
        #pos_dict_file = open(pos_dict_dir, 'w')
        #null_dict_file = open(null_dict_dir, 'w')
        
        with open(tsv_file, 'r') as f:
            header_fields = f.next().strip().split('\t')
            header_fields_copy = header_fields[:13]
            
            #header_fields_copy[13:13] = ['cancer_type', 'peakID', 'peak_type', 'q_value', 'pos_mean', 'pos_count', 'null_mean', 
#                                          'null_count', 'cancer_mean', 'cancer_count', 'normal_mean', 'normal_count',  
#                                          'pos_null_fold_change', 'pos_null_p value', 'pos_normal_fold change', 
#                                          'pos_normal_p_value', 'cancer_normal_fold_change', 'cancer_normal_p_value',
#                                          'plot_hyperlink', 'plot_log_hyperlink', 'hyperlink_antecedent', 'log_hyperlink_antecedent'] 
            cnv_file.write('\t'.join(map(str, header_fields_copy)))
            cnv_file.write('\n')
            reader=csv.reader(f,delimiter='\t')        
            i = 2
            for row in reader:
                pos_expr = []
                null_expr = []
                transcript = row[0]
                peak = transcript_peak[args.toi]
                q_value_peak = peak_info[peak][3]
                #generating count for number of pts with CNV or number of pts w/o CNV
                #don't forget we have the variable "norm_counter" above
                pos_count = len(peak_info[peak][0])
                null_count = len(peak_info[peak][1])
                cancer_count = pos_count + null_count
                
                
                #making list of expression values for each pt with the CNV for this transcript
                
                for x in peak_info[peak][0]:
                    pos_expr.append(float(row[header_fields.index(x)]))
                
                #print to file for plotting 
                #line = [transcript, pos_expr]
                #pos_dict_file.write('\t'.join(map(str, line)))
                #pos_dict_file.write('\n')
                
                #making list of expression values for each pt without the CNV for this transcript
                for x in peak_info[peak][1]:
                    null_expr.append(float(row[header_fields.index(x)]))
                

                
                
#                 #print to file for plotting 
#                 line = [transcript, null_expr]
#                 #null_dict_file.write('\t'.join(map(str, line)))
#                 #null_dict_file.write('\n')
#                 
#                 
#                 #means 
#                 pos_mean = numpy.mean(pos_expr)
#                 null_mean = numpy.mean(null_expr)
#                 normal_mean = norm_dict[transcript]
#                 cancer_mean = (pos_mean*pos_count + null_mean*null_count)/cancer_count
#                 
#                 
#                 #calculate fold changes
#                 pos_null_fc = (pos_mean+0.0000000001)/(null_mean+0.0000000001)
#                 pos_normal_fc = (pos_mean+0.0000000001)/(normal_mean+0.0000000001)
#                 cancer_normal_fc = (cancer_mean+0.0000000001)/(normal_mean+0.0000000001)
#                 
#                 #generate lists for t-test 
#                 cancer_expr = pos_expr + null_expr
#                 normal_expr = norm_expressions_dict[transcript]
#                 
#                 #perform t-tests
#                 t_pos_null,p_pos_null = stats.ttest_ind(pos_expr, null_expr)
#                 t_pos_normal,p_pos_normal = stats.ttest_ind(pos_expr, normal_expr)
#                 t_cancer_normal,p_cancer_normal = stats.ttest_ind(cancer_expr, normal_expr)      
#                 
#                 #create peak-type
#                 peak_type = peak.split('Peak')[0]
                
                #lists of expression values for each category
                #normal --> norm_expressions_dict[transcript]
                #pos --> pos_expr
                #null --> null_expr
                
                #make plots for transcript
                plot_file_dir = os.path.abspath(os.path.join(plot_dir, transcript + '.png'))
                plot_file_log_dir = os.path.abspath(os.path.join(plot_dir, transcript + '_log.png'))
                
                plot_test.plotter(pos_expr, null_expr, norm_expressions_dict[transcript], transcript, plot_file_dir)
                plot_hyperlink = '=HYPERLINK("file://"&AH' + str(i) + '; \"non-log\"'
                plot_test.plotter(pos_expr, null_expr, norm_expressions_dict[transcript], transcript + " log", plot_file_log_dir, True)
                plot_log_hyperlink = '=HYPERLINK("file://"&AI' + str(i) + '; \"log\"'
                
                
#                 row_clip = row[:13]
#                 row_clip[13:13] = [args.tumor_type, peak, peak_type, q_value_peak, pos_mean, pos_count, null_mean, null_count, cancer_mean, cancer_count, 
#                               normal_mean, normal_count, pos_null_fc, p_pos_null, pos_normal_fc, p_pos_normal,
#                               cancer_normal_fc, p_cancer_normal, plot_hyperlink, plot_log_hyperlink, plot_file_dir, plot_file_log_dir ]
                cnv_file.write('\t'.join(map(str, row)))        
                cnv_file.write("\n")
                i+=1
        cnv_file.close()
        #pos_dict_file.close()
        #null_dict_file.close()
    
        

    #run tsv processor for transcript and gene

    
    tsv_processor('gene')
    #tsv_processor('transcript')

    
    p = subprocess.call(['rm', peak_file_dest, peak_file_other_dest, DNA_id_file_dir, lesion_processed_dir, pt_list_dir, normal_ptIDs_dir])    
    p = subprocess.call('rm *.tsv', shell=True)
    

    
    
    
    
    logging.info("Finished")   
    


    return 0



if __name__ == '__main__':

    sys.exit(main())