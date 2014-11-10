'''
Created on Nov 18, 2013

@author: yniknafs
'''

import sys
import argparse
import logging
import collections



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("clinical_file")
        
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    logging.info('Running main script')

    #open clinical file 
    sets_dict = collections.defaultdict(lambda: [])
    null_dict = collections.defaultdict(lambda: [])
    
    clinical_file = open(args.clinical_file, 'r')
    clinical_header = clinical_file.next().strip().split('\t')
    
    universe = []
    
    #generate descriptios for the sets
    basal_dict = collections.defaultdict(lambda: [])
    her2_dict = collections.defaultdict(lambda: [])
    luma_dict = collections.defaultdict(lambda: [])
    lumb_dict = collections.defaultdict(lambda: [])
    normal_dict = collections.defaultdict(lambda: [])
    lumab_dict = collections.defaultdict(lambda: [])
    
    desc_dict = collections.defaultdict(lambda: [])
    desc_dict['PAM_basal'].append('Samples with Basal PAM50 call')
    desc_dict['PAM_her2'].append('Samples with Her2 PAM50RNAseq call')
    desc_dict['PAM_LumA'].append('Samples with LumA PAM50RNAseq call')
    desc_dict['PAM_LumB'].append('Samples with LumB PAM50RNAseq call')
    desc_dict['PAM_Normal'].append('Samples with Normal PAM50 call')
    desc_dict['PAM_LumAB'].append('Samples with LumB PAM50RNAseq call')
    
    
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        sample_type = line[clinical_header.index('sample_type')] 
        if sample_type == "Solid Tissue Normal":
            continue
        universe.append(sampleID)
        #checking stage field 
        stage = line[clinical_header.index('AJCC_Stage_nature2012')]
        
        
        #check PAM50 call
        pam50 = line[clinical_header.index('PAM50Call_RNAseq')]
        # add sample to basal set
        if pam50 == 'Basal':
            sets_dict['PAM_basal'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_basal'].append(sampleID)
        # add sample to HER2 set
        if pam50 == 'Her2':
            sets_dict['PAM_her2'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_her2'].append(sampleID)
        # add sample to LumA set
        if pam50 == 'LumA':
            sets_dict['PAM_LumA'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_LumA'].append(sampleID)
        # add sample to LumB set
        if pam50 == 'LumB':
            sets_dict['PAM_LumB'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_LumB'].append(sampleID)
        # add sample to Normal set
        if pam50 == 'Normal':
            sets_dict['PAM_Normal'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_Normal'].append(sampleID)
        # add sample to lumab set
        if pam50 == 'LumA' or pam50 == 'LumB':
            sets_dict['PAM_LumAB'].append(sampleID)
        elif pam50 != '': 
                null_dict['PAM_LumAB'].append(sampleID)
        
        
    
    smt_file = open('pam50_sets_ab.smt', 'w')
    header_row = ['Set name', 'Set description'] + universe
    smt_file.write('\t'.join(map(str, header_row)))
    smt_file.write('\n')
    for key in sets_dict.iterkeys():
        logging.debug(key)
        desc = desc_dict[key]
        pos = set(sets_dict[key])
        nulls = set(null_dict[key])
        membership = []
        for pt in universe: 
            if pt in pos: 
                membership.append(1)
            elif pt in nulls: 
                membership.append(0)
            else: 
                membership.append('')
        line = [key] + desc + membership
        smt_file.write('\t'.join(map(str,line)))
        smt_file.write('\n')
    smt_file.close()
        

if __name__ == '__main__':
    sys.exit(main())