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
    desc_dict['stage_gte_III'].append('Pts with stage greater than or equal to III')
    desc_dict['stage_IV'].append('Pts with stage of IV')
    desc_dict['metastatic'].append('Pts with metastatic cancer')
    desc_dict['node_positive'].append('Pts with node positive disease')
    desc_dict['node_negative'].append('Pts with node negative disease')
    desc_dict['previous_nodular_hyperplasia'].append('Pts with nodular hyperplasia previous to cancer diagnosis')
    desc_dict['previous_lymphocytic_thyroiditis'].append('Pts with lymphocytic_thyroiditis previous to cancer diagnosis')
    desc_dict['unifocal'].append('Pts with unifocal disease')
    desc_dict['multifocal'].append('Pts with multifocal disease')
    
    

            
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)
        #checking age at diagnosis field 
                #check Days to date of death field
        dtd = line[clinical_header.index('_OS')]
        dead = line[clinical_header.index('_OS_IND')]
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
        
        #add to stage >=3 set
        if stage == 'Stage III' or \
            stage == 'Stage IVA' or \
            stage == 'Stage IVB' or \
            stage == 'Stage IVC':
            sets_dict['stage_gte_III'].append(sampleID)
        elif stage: 
            null_dict['stage_gte_III'].append(sampleID)
            
        
        #add to stage =4 set
        if  stage == ('Stage IV' or 'Stage IVA' or 'Stage IVC'):        
            sets_dict['stage_IV'].append(sampleID)
        elif stage: 
            null_dict['stage_IV'].append(sampleID)    
        
        
        met = line[clinical_header.index('pathologic_M')]
        #add to metastatic set
        if met == 'M1':
            sets_dict['metastatic'].append(sampleID)
        elif met: 
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
        
        #check node status
        disease = line[clinical_header.index('patient_personal_medical_history_thyroid_gland_disorder_name')]
        if disease == 'Nodular Hyperplasia':
            sets_dict['previous_nodular_hyperplasia'].append(sampleID)
        elif disease: 
            null_dict['previous_nodular_hyperplasia'].append(sampleID)
        if disease == 'Lymphocytic Thyroiditis':
            sets_dict['previous_lymphocytic_thyroiditis'].append(sampleID)
        elif disease: 
            null_dict['previous_lymphocytic_thyroiditis'].append(sampleID)

        #check node status
        focus = line[clinical_header.index('primary_neoplasm_focus_type')]
        if focus == 'Multifocal':
            sets_dict['multifocal'].append(sampleID)
        elif focus: 
            null_dict['multifocal'].append(sampleID)
        if focus == 'Unifocal':
            sets_dict['unifocal'].append(sampleID)
        elif focus: 
            null_dict['unifocal'].append(sampleID)
        
        
    ca_type = 'THCA' 
    
    smt_file = open(ca_type + '_clinical_sets.smt', 'w')
    header_row = ['Set name', 'Set description'] + universe
    smt_file.write('\t'.join(map(str, header_row)))
    smt_file.write('\n')
    for key in sets_dict.iterkeys():
        print key
        print desc_dict[key]
        
        desc = ca_type + ' ' +desc_dict[key][0]
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
        name = ca_type + ' ' + key.replace('_', ' ')
        line = [name] + [desc] + membership
        smt_file.write('\t'.join(map(str,line)))
        smt_file.write('\n')
    smt_file.close()
        

if __name__ == '__main__':
    sys.exit(main())