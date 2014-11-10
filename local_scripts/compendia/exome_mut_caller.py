import argparse
import logging
import os
import sys
from sets import Set
import collections

'''

'''



def read_lines(filename, header = True):
    lines = []
    i = 0
    for line in open(filename):
        if header == True and i == 0:            
            i += 1
            header_line = line
            continue 
        if not line:
            continue
        if line.startswith("#"):
            continue
        line = line.strip()
        if not line:
            continue
        lines.append(line)
        
    if header == True: 
        return lines, header_line    
    else:
        return lines

def main():
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("phenos")
    parser.add_argument("mut_file")
    parser.add_argument("meta_file")
    parser.add_argument("output")
    parser.add_argument("-c", dest = "compendia_ID_file",
                        default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    args = parser.parse_args()
    
    
    
    #read file that contains two columns: compendiaID, TCGAsampleID 
    compendia_IDs, header = read_lines(args.compendia_ID_file)
    #make a dictionary to convert TCGAptID to the compendia sampleIDs
    id_conv = collections.defaultdict(lambda: 'NOT_IN_COMPENDIA')
    for line in compendia_IDs:
        line_items = line.split('\t')
        if len(line_items)>1:
            libID = line_items[0]
            TCGAID = line_items[1]
            #ID_file contains non-TCGA samples
            if TCGAID.startswith("TCGA"):
                    #the TCGAptID has more numbers than the clinicalMatrix file
                    #so the ptID is clipped to match length
                    TCGA_split = TCGAID.split('-')
                    #remove the letter from the end of the last hyphen separated item
                    TCGA_split3 = TCGA_split[3][:2]
                    #join the first 4 items separated by hyphen
                    
                    TCGA_newID = '-'.join([TCGA_split[0], TCGA_split[1], TCGA_split[2], TCGA_split3])
                    id_conv[TCGA_newID] = libID
    
    #make dictionary matching sampleID to TCGA cancer type
    ca_type_dict = collections.defaultdict(lambda: 'NOT_IN_COMPENDIA')
    file = open(args.phenos, 'r')
    header = file.next().strip().split('\t')
    for line in file: 
        line = line.strip().split('\t')
        libID = line[header.index('library_id')]
        TCGA_ca = line[header.index('tcga_disease_type')]
        ca_type_dict[libID] = TCGA_ca
    
    gene_len_dict = collections.defaultdict(lambda: 'NOT_IN_COMPENDIA')
    file = open(args.meta_file, 'r')
    header = file.next().strip().split('\t')
    for line in file: 
        line = line.strip().split('\t')
        gene = line[header.index('gene_name')]
        length = line[header.index('transcript_length')]
        gene_len_dict[gene] = int(length)
    
    file = open(args.mut_file, 'r')
    header = file.next().strip().split('\t')
    #convert to compendia lib IDs and count pts in each cancer type
    ca_count = collections.defaultdict(lambda: [])
    for x in xrange(1, len(header)): 
        libID_conv = id_conv[header[x]]
        header[x] = libID_conv
        ca_type = ca_type_dict[libID_conv]
        ca_count[ca_type].append(libID_conv)
        ca_count['pancan'].append(libID_conv)
#     for key in ca_count: 
#         print key, len(ca_count[key])
    
    fileo1 = open('muts.tsv', 'w')
    fileo2 = open('muts_norm_len.tsv', 'w')
    fileo3 = open('muts_norm_len_type.tsv', 'w')
    keys = []
    headers = []
    for key in ca_count.iterkeys(): 
        headers.append(key + ('(n=%d)' % len(ca_count[key])))
        keys.append(key)
        
    headero = ['gene', 'length'] + headers
    fileo1.write('\t'.join(headero)+'\n')
    fileo2.write('\t'.join(headero)+'\n')
    fileo3.write('\t'.join(headero)+'\n')
    
    tot = 0 
    for line in file: 
        tot +=1
    
    file = open(args.mut_file, 'r')
    file.next()
    i = 0
    logging.debug("Counting mutations")
    for line in file:
        i+=1
        if (i%2500) == 0: 
            logging.debug('Mutations counted for %d/%d genes' % (i, tot))
        count_dict = collections.defaultdict(lambda: [])
        line = line.strip().split('\t')
        gene = line[header.index('Sample')]
        length = gene_len_dict[gene]
        if length == "NOT_IN_COMPENDIA":
            print gene
            continue
        if length < 200:
            print gene
            continue
        for x in xrange(1, len(header)):
            ca_type = ca_type_dict[header[x]]
            count_dict[ca_type].append(float(line[x]))
            count_dict['pancan'].append(float(line[x]))
        lineo1 = [gene, length]
        lineo2 = [gene, length]        
        lineo3 = [gene, length]
        for key in keys: 
            num_muts = sum(count_dict[key])
            num_muts_norm = num_muts/length
            num_muts_norm_type = num_muts/length/(len(ca_count[key])) 
            lineo2.append(num_muts_norm)
            lineo1.append(num_muts)
            lineo3.append(num_muts_norm_type)
        fileo1.write('\t'.join(map(str,lineo1)) + '\n')
        fileo2.write('\t'.join(map(str,lineo2)) + '\n')
        fileo3.write('\t'.join(map(str,lineo3)) + '\n')
        
     
    
    
    
    
if __name__ == '__main__':
    sys.exit(main())