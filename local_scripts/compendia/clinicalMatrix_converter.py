import argparse
import logging
import os
import sys
from sets import Set
import collections

'''
Takes the clinical matrix file downloaded from the 
UCSC cancer genome browser and converts the TCGA
ptIDs into libraryIDs that the compendia assembly 
uses. 
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
    parser.add_argument("clinical_matrix")
    parser.add_argument("output")
    parser.add_argument("-c", dest = "compendia_ID_file",
                        default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    parser.add_argument("--cnv", dest = "cnv",
                        default = None)
    args = parser.parse_args()
    
    
    if args.cnv != None:
        cnv_IDs = []
        cnv_file = open(args.cnv, 'r')
        for line in cnv_file:
            line = line.strip()
            cnv_IDs.append(line)
        cnv_IDs = set(cnv_IDs)
    
    
    #read file that contains two columns: compendiaID, TCGAsampleID 
    compendia_IDs, header = read_lines(args.compendia_ID_file)
    #make a dictionary to convert TCGAptID to the compendia sampleIDs
    id_conv = collections.defaultdict(lambda: [])
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
    
    
    
    matrix_lines, header_matrix = read_lines(args.clinical_matrix) 
    out_file = open(args.output + '.txt', 'w')
    out_file.write(header_matrix)
    for line in matrix_lines:
        line_items = line.split('\t') 
        TCGA_id = line_items[0]
        #check to see if this ptID is in the compendia
        converted_id = id_conv[TCGA_id] 
        if converted_id != []:
            if args.cnv != None:
                if converted_id in cnv_IDs: 
                    line_items[0] = id_conv[TCGA_id]
                    out_file.write('\t'.join(map(str, line_items)))
                    out_file.write('\n')
            else: 
                line_items[0] = id_conv[TCGA_id]
                out_file.write('\t'.join(map(str, line_items)))
                out_file.write('\n')
    out_file.close()
    
    return 0


#def clinicalMatrix_reader(file, header = True):

    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    sys.exit(main())