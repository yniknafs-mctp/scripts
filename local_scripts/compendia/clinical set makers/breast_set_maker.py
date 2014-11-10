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
    desc_dict = collections.defaultdict(lambda: [])
    desc_dict['stage_gte_II'].append('Samples with stage greater than or equal to II')
    desc_dict['stage_gte_III'].append('Samples with stage greater than or equal to III')
    desc_dict['metastatic'].append('Samples from pts with metastasis')
    desc_dict['dead_3yrs'].append('Pts who died within three years of diagnosis')
    desc_dict['dead_5yrs'].append('Pts who died within five years of diagnosis')    
    desc_dict['dead_10yrs'].append('Pts who died within ten years of diagnosis')
    desc_dict['ER_positive'].append('Samples that are ER positive')
    desc_dict['ER_negative'].append('Samples that are ER negative')
    desc_dict['HER2_positive'].append('Samples that are HER2 positive')
    desc_dict['HER2_negative'].append('Samples that are HER2 negative')
    desc_dict['PR_positive'].append('Samples that are PR positive')
    desc_dict['PR_negative'].append('Samples that are PR negative')
    desc_dict['triple_negative'].append('Samples that are triple negative')
    desc_dict['node_positive'].append('Samples that are node positive')
    desc_dict['node_negative'].append('Samples that are node negative')
    desc_dict['PAM_basal'].append('Samples with Basal PAM50 call')
    desc_dict['PAM_her2'].append('Samples with Her2 PAM50RNAseq call')
    desc_dict['PAM_LumA'].append('Samples with LumA PAM50RNAseq call')
    desc_dict['PAM_LumB'].append('Samples with LumB PAM50RNAseq call')
    desc_dict['PAM_Normal'].append('Samples with Normal PAM50 call')
    
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)
        #checking stage field 
        stage = line[clinical_header.index('AJCC_Stage_nature2012')]
        
        #add to stage >=2 set
        if stage == 'Stage II' or \
            stage == 'Stage IIA' or \
            stage == 'Stage IIB' or \
            stage == 'Stage III' or \
            stage == 'Stage IIIA' or \
            stage == 'Stage IIIB' or \
            stage == 'Stage IIIC' or \
            stage == 'Stage IV':
            sets_dict['stage_gte_II'].append(sampleID)
        elif stage != '': 
            null_dict['stage_gte_II'].append(sampleID)
            
        
        #add to stage >=3 set
        if  stage == 'Stage III' or \
            stage == 'Stage IIIA' or \
            stage == 'Stage IIIB' or \
            stage == 'Stage IIIC' or \
            stage == 'Stage IV':        
            sets_dict['stage_gte_III'].append(sampleID)
        elif stage != '': 
            null_dict['stage_gte_III'].append(sampleID)    
            
        #add to metastatic set
        if stage == 'Stage IV':
            sets_dict['metastatic'].append(sampleID)
        elif stage != '': 
            null_dict['metastatic'].append(sampleID)
            
        #check Days to date of death field
        dtd = line[clinical_header.index('Days_to_date_of_Death_nature2012')]
        dtfu = line[clinical_header.index('Days_to_Date_of_Last_Contact_nature2012')]
        if dtd: 
            #add pts to dead at 3yrs
            if float(dtd) < 1095: 
                sets_dict['dead_3yrs'].append(sampleID)
            else: 
                null_dict['dead_3yrs'].append(sampleID)
            #add pts to dead at 5yrs
            if float(dtd) < 1825:
                sets_dict['dead_5yrs'].append(sampleID)
            else:
                null_dict['dead_5yrs'].append(sampleID)
            #add pts to dead at 10yrs
            if float(dtd) < 3650:
                sets_dict['dead_10yrs'].append(sampleID)
            else:
                null_dict['dead_10yrs'].append(sampleID)
        else: 
            if dtfu: 
                if float(dtfu) >= 1095:
                    null_dict['dead_3yrs'].append(sampleID)
                if float(dtfu) >= 1825:
                    null_dict['dead_5yrs'].append(sampleID)
                if float(dtfu) >= 3650:
                    null_dict['dead_10yrs'].append(sampleID)
        #check ER status
        ER = line[clinical_header.index('ER_Status_nature2012')]
        #add pts to ER positive/negative set
        if ER == "Positive": 
            sets_dict['ER_positive'].append(sampleID)
        elif ER != '': 
                null_dict['ER_positive'].append(sampleID)
        if ER == "Negative":
            sets_dict['ER_negative'].append(sampleID)
        elif ER != '': 
                null_dict['ER_negative'].append(sampleID)
        
        
        #check HER2 status
        HER2 = line[clinical_header.index('HER2_Final_Status_nature2012')]
        # add pts to HER2 positive/negative set
        if HER2 == "Positive": 
            sets_dict['HER2_positive'].append(sampleID)
        elif HER2 != '': 
                null_dict['HER2_positive'].append(sampleID)
        if HER2 == "Negative":
            sets_dict['HER2_negative'].append(sampleID)
        elif HER2 != '': 
                null_dict['HER2_negative'].append(sampleID)
        
        #check PR status
        PR = line[clinical_header.index('PR_Status_nature2012')]
        # add pts to HER2 positive/negative set
        if PR == "Positive": 
            sets_dict['PR_positive'].append(sampleID)
        elif PR != '': 
                null_dict['PR_positive'].append(sampleID)
        if PR == "Negative":
            sets_dict['PR_negative'].append(sampleID)
        elif PR != '': 
                null_dict['PR_negative'].append(sampleID)
            
        #make triple negative set
        if ER == "Negative" and HER2 == "Negative" and PR == "Negative":
            sets_dict['triple_negative'].append(sampleID)
        elif PR != '' or ER != '' or HER2 != '': 
                null_dict['triple_negative'].append(sampleID)
        
        #check node status
        node = line[clinical_header.index('Node_Coded_nature2012')]
        # add pts to HER2 positive/negative set
        if node == "Positive": 
            sets_dict['node_positive'].append(sampleID)
        elif node != '': 
                null_dict['node_positive'].append(sampleID)
        if node == "Negative":
            sets_dict['node_negative'].append(sampleID)
        elif node != '': 
                null_dict['node_negative'].append(sampleID)
            
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
    
    smt_file = open('breast_sets.smt', 'w')
    header_row = ['Set name', 'Set description'] + universe
    smt_file.write('\t'.join(map(str, header_row)))
    smt_file.write('\n')
    for key in sets_dict.iterkeys():
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