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
    desc_dict['dead_3yrs'].append('Pts who were dead at 3 years after diagnosis')
    desc_dict['dead_5yrs'].append('Pts who were dead at 5 years after diagnosis')
    desc_dict['dead_10yrs'].append('Pts who were dead at 10 years after diagnosis')
    desc_dict['male'].append('male patients')
    desc_dict['female'].append('female patients')
    desc_dict['stage_gte_II'].append('Pts with stage greater than or equal to II')
    desc_dict['stage_gte_III'].append('Pts with stage greater than or equal to III')
    desc_dict['breslow_gte_1.5'].append('Breslow depth greater than or equal to 1.5mm (Breslow stage III and greater)')
    desc_dict['breslow_gte_3'].append('Breslow depth greater than or equal to 3.0mm (Breslow stage V)')
    desc_dict['mitotic_rate_lte_4'].append('Mitotic rate less than 4 per square mm')
    desc_dict['mitotic_rate_gte_4'].append('Mitotic rate greater than 4 per square mm')
    desc_dict['Clark_level_gte_III'].append('Clark level greater than or equal to III')
    desc_dict['Clark_level_gte_IV'].append('Clark level greater than or equal to III')
    desc_dict['Clark_level_V'].append('Clark level of V')
    desc_dict['ulceration_positive'].append('Pts with ulcerated lesions')
    desc_dict['ulceration_negative'].append('Pts with ulcerated lesions')
    desc_dict['metastatic'].append('Pts with metastatic cancer')
    desc_dict['recurrence'].append('Pts with locoregional recurrence')
    desc_dict['node_positive'].append('Pts with node positive disease')
    desc_dict['node_negative'].append('Pts with node negative disease')
    desc_dict['interferon'].append('Pts who were treated with interferon prior to resection')
    
    

            
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)
        #checking age at diagnosis field 
                #check Days to date of death field
        dtd = line[clinical_header.index('_OS')]
        dead = line[clinical_header.index('_EVENT')]
        if dtd: 
            #add pts to dead at 3yrs
            if float(dtd) < 1095 and dead == '1': 
                sets_dict['dead_3yrs'].append(sampleID)
            elif dtd: 
                null_dict['dead_3yrs'].append(sampleID)
            #add pts to dead at 5yrs
            if float(dtd) < 1825 and dead == '1':
                sets_dict['dead_5yrs'].append(sampleID)
            elif dtd: 
                null_dict['dead_5yrs'].append(sampleID)
            #add pts to dead at 10yrs
            if float(dtd) < 3650 and dead == '1':
                sets_dict['dead_10yrs'].append(sampleID)
            elif dtd:
                null_dict['dead_10yrs'].append(sampleID)
                
        stage = line[clinical_header.index('pathologic_stage')]
        
        #add to stage >=2 set
        if stage == 'Stage II' or \
            stage == 'Stage IIA' or \
            stage == 'Stage IIB' or \
            stage == 'Stage IIC' or \
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
        
        met = line[clinical_header.index('Distant Metastasis')]
        #add to metastatic set
        if met == 'Yes':
            sets_dict['metastatic'].append(sampleID)
        elif met != '': 
            null_dict['metastatic'].append(sampleID)
        
        
        #check gender
        gender = line[clinical_header.index('gender')]
        #add pts to male
        if gender == 'MALE': 
            sets_dict['male'].append(sampleID)
        elif gender: 
            null_dict['male'].append(sampleID)
        #add pts to female
        if gender == 'FEMALE': 
            sets_dict['female'].append(sampleID)
        elif gender: 
            null_dict['female'].append(sampleID)
            
        #check breslow
        breslow = line[clinical_header.index('breslow_depth_value')]
        if breslow:
            breslow = float(breslow)
            if breslow > 1.5: 
                sets_dict['breslow_gte_1.5'].append(sampleID)
            elif breslow: 
                null_dict['breslow_gte_1.5'].append(sampleID)
            if breslow > 3: 
                sets_dict['breslow_gte_3'].append(sampleID)
            elif breslow: 
                null_dict['breslow_gte_3'].append(sampleID)
        
        #check mr
        mr = line[clinical_header.index('malignant_neoplasm_mitotic_count_rate')]
        if mr:
            mr = float(mr)
            if mr < 4: 
                sets_dict['mitotic_rate_lte_4'].append(sampleID)
            elif mr: 
                null_dict['mitotic_rate_lte_4'].append(sampleID)
            if mr > 4: 
                sets_dict['mitotic_rate_gte_4'].append(sampleID)
            elif mr: 
                null_dict['mitotic_rate_gte_4'].append(sampleID)
        
        
        #check cl
        cl = line[clinical_header.index('melanoma_clark_level_value')]
        if cl == 'III' or \
                cl == 'IV' or \
                cl == 'V': 
            sets_dict['Clark_level_gte_III'].append(sampleID)
        elif cl: 
            null_dict['Clark_level_gte_III'].append(sampleID)
        if cl == 'IV' or \
                cl == 'V': 
            sets_dict['Clark_level_gte_IV'].append(sampleID)
        elif cl: 
            null_dict['Clark_level_gte_IV'].append(sampleID)
        if cl == 'V': 
            sets_dict['Clark_level_V'].append(sampleID)
        elif cl: 
            null_dict['Clark_level_V'].append(sampleID)
        
        #check ulceration
        ulcer = line[clinical_header.index('melanoma_ulceration_indicator')]
        if ulcer == 'YES':
            sets_dict['ulceration_positive'].append(sampleID)
        elif ulcer: 
            null_dict['ulceration_positive'].append(sampleID)
        if ulcer == 'NO':
            sets_dict['ulceration_negative'].append(sampleID)
        elif ulcer: 
            null_dict['ulceration_negative'].append(sampleID)
        
        #check recurrence
        recurrence = line[clinical_header.index('Locoregional Recurrence')]
        if recurrence == 'Yes':
            sets_dict['recurrence'].append(sampleID)
        elif recurrence: 
            null_dict['recurrence'].append(sampleID)
        
        #check node status
        node = line[clinical_header.index('pathologic_N')]
        if node != ('N0' or 'NX' or ''):
            sets_dict['node_positive'].append(sampleID)
        elif node: 
            null_dict['node_positive'].append(sampleID)
        if node == ('N0' or 'NX'):
            sets_dict['node_negative'].append(sampleID)
        elif node: 
            null_dict['node_negative'].append(sampleID)
        
        
    
    smt_file = open('SKCM_clinical_sets.smt', 'w')
    header_row = ['Set name', 'Set description'] + universe
    smt_file.write('\t'.join(map(str, header_row)))
    smt_file.write('\n')
    for key in sets_dict.iterkeys():
        desc = 'SKCM ' +desc_dict[key][0]
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
        name = 'SKCM ' + key.replace('_', ' ')
        line = [name] + [desc] + membership
        smt_file.write('\t'.join(map(str,line)))
        smt_file.write('\n')
    smt_file.close()
        

if __name__ == '__main__':
    sys.exit(main())