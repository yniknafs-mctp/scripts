import argparse
import logging
import os
import sys
from sets import Set
import collections

'''
Takes any file with columns representing TCGA pt_IDs and 
converts the TCGA IDs to RNAseq IDs

currently only supports columns conversion
'''



def main():
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("file_to_convert")
    parser.add_argument("--col", dest = 'col',
                        default = False)
    parser.add_argument("--row", dest = 'row',
                        default = False)
    parser.add_argument("-c", dest = "compendia_ID_file",
                        default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    args = parser.parse_args()
    
        
    #read file that contains two columns: compendiaID, TCGAsampleID 
    compendia_ID_file = open(args.compendia_ID_file)
    compendia_ID_file.next()
    id_conv = collections.defaultdict(lambda: [])
    for line in compendia_ID_file: 
        line = line.strip().split('\t')
        if len(line)>1:
            libID = line[0]
            TCGAID = line[1]
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
    
    
    #read file to convert 
    file_to_convert = open(args.file_to_convert)
    header = file_to_convert.next().strip().split("\t")
    if args.col:
        for x in xrange(int(args.col), len(header)):
            conversion = id_conv[header[x]]
            if conversion != []:
                header[x] = conversion
            else: 
                header[x] = 'NOT_IN_COMPENDIA' 
        print '\t'.join(header)
        for line in file_to_convert: 
            print line.strip()

    
    
    return 0

if __name__ == '__main__':
    sys.exit(main())