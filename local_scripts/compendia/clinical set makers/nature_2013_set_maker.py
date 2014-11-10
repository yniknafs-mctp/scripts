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

LAST_META_COL = 12
FRAC_CUTOFF = 0.05
NUM_SAMPLES_CUTOFF = 8

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
    parser.add_argument("mut_matrix_file")
    parser.add_argument("phenos")
    parser.add_argument("-c", dest = "compendia_ID_file",
                        default = '/home/yniknafs/misc_files/CNV_ptID_long_to_libID.txt')
    args = parser.parse_args()
    
    
    #read file that contains two columns: compendiaID, TCGAsampleID 
    compendia_IDs, header = read_lines(args.compendia_ID_file)
    #make a dictionary to convert TCGAptID to the compendia sampleIDs
    id_conv = collections.defaultdict(lambda: "NOT_IN_COMPENDIA")
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
                    
                    TCGA_newID = '-'.join([TCGA_split[0], TCGA_split[1], TCGA_split[2]])
                    id_conv[TCGA_newID] = libID
    
    ca_type_dict = collections.defaultdict(lambda: "NOT_IN_COMPENDIA")
    file = open(args.phenos, 'r')
    header = file.next().strip().split('\t')
    for line in file: 
        line = line.strip().split('\t')
        libID = line[header.index('library_id')]
        TCGA_ca = line[header.index('tcga_disease_type')]
        if TCGA_ca == ("COAD" or "READ"):
            TCGA_ca = "COADREAD"
        ca_type_dict[libID] = TCGA_ca
    
    #read in matrix file 
    matrix_lines, header_matrix = read_lines(args.mut_matrix_file)
    header_matrix = header_matrix.strip().split('\t')
    genes = header_matrix[LAST_META_COL:]
    #make new lists that convert the gene name and description field into readable format
    def pancan_namer(gene):
        prefix = 'Pancan'
        return prefix + " " + gene + ' mutation'
    def pancan_describer(gene):
        prefix = 'Pancan'
        return "%s patients with a somatic mutation in %s" % (prefix, gene)
    names = map(pancan_namer, genes)
    desc = map(pancan_describer, genes)
    
    
    fileo = open('mutation_sets.smx', 'w')
    #write file first and second line (set name, set desc)
    fileo.write('\t'.join(['Set name'] + names) + '\n')
    fileo.write('\t'.join(['Set description'] + desc) + '\n')
    
    #create dictionaries to count mutation data for each cancer type
    count_tot = collections.defaultdict(lambda: 0)
    count_muts = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    hits_muts = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    
    #make list of all samples to be used to print cancer type sets
    universe = []
    for line in matrix_lines:
        line = line.strip().split('\t')
        #split sampleID to be paired to the compendia key
        sample = line[0]
        sample_split = sample.split('-')
        sample_strip = '-'.join([sample_split[0], sample_split[1], sample_split[2]])
        compID = id_conv[sample_strip]
        #keep list of all sampleIDs that are in compendia
        if compID != "NOT_IN_COMPENDIA":
            universe.append(compID)
        #get cancer type for sample
        ca_type = ca_type_dict[compID]
        #skip samples not in the compendia
        if compID == "NOT_IN_COMPENDIA" or ca_type == "NOT_IN_COMPENDIA" or ca_type == 'LAML':
            continue
        #add to the count_tot dict to get the total number of samples for each type
        count_tot[ca_type]+=1
        #loop through all the genes and add to the lists for each cancer type
        line_genes = line[LAST_META_COL:]
        for i in xrange(len(genes)): 
            gene = genes[i]
            mut_call = int(line_genes[i])
            if mut_call == 1:
                hits_muts[ca_type][gene].append(compID)
                count_muts[ca_type][gene]+=1
                
        #Print to pancan set file
        calls = line[LAST_META_COL:]
        lineo = [compID] + calls
        fileo.write('\t'.join(lineo)+'\n')
    
    #loop through the cancer type dictionaries to make sets that meet cutoffs
    fileo2 = open('mutation_type_sets.smt', 'w')
    header = ['Sample name', 'Sample description'] + universe
    fileo2.write('\t'.join(header) + '\n')
    for type_key in hits_muts.iterkeys():
        for gene_key in hits_muts[type_key].iterkeys():
            count = float(count_muts[type_key][gene_key]) 
            tot = float(count_tot[type_key])
            frac = count/tot
            if frac > FRAC_CUTOFF and count > NUM_SAMPLES_CUTOFF: 
                set_name = "%s %s mutation" % (type_key, gene_key)
                set_desc = "%s patients with a somatic mutation in %s" % (type_key, gene_key)
                #loop through universe array to determine membership 
                membership = []
                hits = set(hits_muts[type_key][gene_key])
                for sample in universe: 
                    if sample in hits: 
                        membership.append(1)
                    elif ca_type_dict[sample] == type_key: 
                        membership.append(0)
                    else: 
                        membership.append('')
                lineo = [set_name, set_desc] + membership
                fileo2.write('\t'.join(map(str, lineo)) + '\n')
                    
                
    fileo.close()
    fileo2.close()
    return 0


#def clinicalMatrix_reader(file, header = True):

    
    
    
    
    
    
    
    
    
    
if __name__ == '__main__':
    sys.exit(main())