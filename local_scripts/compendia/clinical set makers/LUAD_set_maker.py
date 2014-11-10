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
    desc_dict['BRAF_point_mutation'].append('patients with a BRAF point mutation')
    desc_dict['CBL_point_mutation'].append('patients with a CBL point mutation')
    desc_dict['CTNNB1_point_mutation'].append('patients with a CTNNB1 point mutation')
    desc_dict['KRAS_EGFR_ALK_point_mutations'].append('patients with multiple mutations in KRAS, EGFR, and ALK')
    desc_dict['KRAS_EGFR_ALK_RET_ROS1_BRAF_ERBB2_HRAS_NRAS_AKT1_MAP2_point_mutations'].append('patients with multiple mutations in KRAS, EGFR, ALK, RET, ROS1, BRAF, \
                                                                                                    ERBB2, HRAS, NRAS, AKT1, and MAP2')
    desc_dict['EGFR_point_mutation'].append('patients with an EGFR point mutation')
    desc_dict['ERBB2_point_mutation'].append('patients with an ERBB2 point mutation')
    desc_dict['ERBB4_point_mutation'].append('patients with an ERBB4 point mutation')
    desc_dict['HRAS_point_mutation'].append('patients with an HRAS point mutation')
    desc_dict['KRAS_point_mutation'].append('patients with a KRAS point mutation')
    desc_dict['NRAS_point_mutation'].append('patients with an NRAS point mutation')
    desc_dict['PIK3CA_point_mutation'].append('patients with a PIK3CA point mutation')
    desc_dict['STK11_point_mutation'].append('patients with a STK11 point mutation')
    desc_dict['MAP2K1_point_mutation'].append('patients with a MAP2K1 point mutation')
    desc_dict['PTPN11_point_mutation'].append('patients with a PTPN11 point mutation')
    desc_dict['MET_point_mutation'].append('patients with a MET point mutation')
    desc_dict['expression_subtype_bronchoid'].append('patients with bronchoid expression subtype')
    desc_dict['expression_subtype_magnoid'].append('patients with magnoid expression subtype')
    desc_dict['expression_subtype_squamoid'].append('patients with squamoid expression subtype')
    desc_dict['dead_3yrs'].append('Pts who were dead at 3 years after diagnosis')
    desc_dict['dead_5yrs'].append('Pts who were dead at 5 years after diagnosis')
    desc_dict['dead_1yr'].append('Pts who were dead at 1 year after diagnosis')
    desc_dict['male'].append('male patients')
    desc_dict['female'].append('female patients')
    desc_dict['node_positive'].append('Pts with node positive disease')
    desc_dict['node_negative'].append('Pts with node negative disease')
    desc_dict['stage_gte_II'].append('Pts with stage greater than or equal to II')
    desc_dict['stage_gte_III'].append('Pts with stage greater than or equal to III')
    desc_dict['progressive_disease'].append('Pts with progressive disease')
    desc_dict['disease_in_remission'].append('Pts whose disease is in remission')
    desc_dict['reformed_smoker_lte_15yrs'].append('Pts who quit smoking less than 15 years ago')
    desc_dict['reformed_smoker_gt_15yrs'].append('Pts who quit smoking more than 15 years ago')
    desc_dict['current_smoker'].append('Pts who currently smoke')
    desc_dict['never_smoker'].append('Pts who are lifelong nonsmokers')
    desc_dict['metastatic'].append('Pts with metastatic cancer')
    desc_dict['radiation_therapy'].append('Pts who received radiation therapy prior to resection')
    
       
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)

        braf = line[clinical_header.index('BRAF')]
        if braf != ('none') and braf != '' :
            print braf
            sets_dict['BRAF_point_mutation'].append(sampleID)
        elif braf: 
            null_dict['BRAF_point_mutation'].append(sampleID)

        cbl = line[clinical_header.index('CBL')]
        if cbl != ('none') and cbl != '':
            sets_dict['CBL_point_mutation'].append(sampleID)
        elif cbl: 
            null_dict['CBL_point_mutation'].append(sampleID)

        ctnnb1 = line[clinical_header.index('CTNNB1')]
        if ctnnb1 != ('none') and ctnnb1 != '' :
            sets_dict['CTNNB1_point_mutation'].append(sampleID)
        elif ctnnb1: 
            null_dict['CTNNB1_point_mutation'].append(sampleID)

        KRAS_EGFR_ALK = line[clinical_header.index('Canonical_mut_in_KRAS_EGFR_ALK')]
        if KRAS_EGFR_ALK == 'Y':
            sets_dict['KRAS_EGFR_ALK_point_mutations'].append(sampleID)
        elif KRAS_EGFR_ALK: 
            null_dict['KRAS_EGFR_ALK_point_mutations'].append(sampleID)

        mut_combo = line[clinical_header.index('Cnncl_mt_n_KRAS_EGFR_ALK_RET_ROS1_BRAF_ERBB2_HRAS_NRAS_AKT1_MAP2')]
        if mut_combo == 'Y':
            sets_dict['KRAS_EGFR_ALK_RET_ROS1_BRAF_ERBB2_HRAS_NRAS_AKT1_MAP2_point_mutations'].append(sampleID)
        elif mut_combo: 
            null_dict['KRAS_EGFR_ALK_RET_ROS1_BRAF_ERBB2_HRAS_NRAS_AKT1_MAP2_point_mutations'].append(sampleID)

        EGFR = line[clinical_header.index('EGFR')]
        if EGFR != ('none') and EGFR != '':
            sets_dict['EGFR_point_mutation'].append(sampleID)
        elif EGFR: 
            null_dict['EGFR_point_mutation'].append(sampleID)
        
        ERBB2 = line[clinical_header.index('ERBB2')]
        if ERBB2 != ('none') and ERBB2 != '' :
            sets_dict['ERBB2_point_mutation'].append(sampleID)
        elif ERBB2: 
            null_dict['ERBB2_point_mutation'].append(sampleID)

        ERBB4 = line[clinical_header.index('ERBB4')]
        if ERBB4 != ('none') and ERBB4 !='':
            sets_dict['ERBB4_point_mutation'].append(sampleID)
        elif ERBB4: 
            null_dict['ERBB4_point_mutation'].append(sampleID)

        HRAS = line[clinical_header.index('HRAS')]
        if HRAS != ('none') and HRAS !='':
            sets_dict['HRAS_point_mutation'].append(sampleID)
        elif HRAS: 
            null_dict['HRAS_point_mutation'].append(sampleID)

        KRAS = line[clinical_header.index('KRAS')]
        if KRAS != ('none') and KRAS !='' :
            sets_dict['KRAS_point_mutation'].append(sampleID)
        elif KRAS: 
            null_dict['KRAS_point_mutation'].append(sampleID)        
        
        NRAS = line[clinical_header.index('NRAS')]
        if NRAS != ('none') and NRAS!='':
            sets_dict['NRAS_point_mutation'].append(sampleID)
        elif NRAS: 
            null_dict['NRAS_point_mutation'].append(sampleID)
        
        PIK3CA = line[clinical_header.index('PIK3CA')]
        if PIK3CA != ('none') and PIK3CA!='':
            sets_dict['PIK3CA_point_mutation'].append(sampleID)
        elif PIK3CA: 
            null_dict['PIK3CA_point_mutation'].append(sampleID)
        
        PTPN11 = line[clinical_header.index('PTPN11')]
        if PTPN11 != ('none') and PTPN11!='':
            sets_dict['PTPN11_point_mutation'].append(sampleID)
        elif PTPN11: 
            null_dict['PTPN11_point_mutation'].append(sampleID)
        
        STK11 = line[clinical_header.index('STK11')]
        if STK11 != ('none') and STK11!='':
            sets_dict['STK11_point_mutation'].append(sampleID)
        elif STK11: 
            null_dict['STK11_point_mutation'].append(sampleID)
            
        MAP2K1 = line[clinical_header.index('MAP2K1')]
        if MAP2K1 != ('none') and (MAP2K1 !='') :
            print MAP2K1
            sets_dict['MAP2K1_point_mutation'].append(sampleID)
        elif MAP2K1: 
            null_dict['MAP2K1_point_mutation'].append(sampleID)
            
        MET = line[clinical_header.index('MET')]
        if MET != ('none') and MET !='' :
            sets_dict['MET_point_mutation'].append(sampleID)
        elif MET: 
            null_dict['MET_point_mutation'].append(sampleID)
            
        expr = line[clinical_header.index('Expression_Subtype')]
        if expr == 'Bronchoid':
            sets_dict['expression_subtype_bronchoid'].append(sampleID)
        elif expr: 
            null_dict['expression_subtype_bronchoid'].append(sampleID)
        
        expr = line[clinical_header.index('Expression_Subtype')]
        if expr == 'Magnoid':
            sets_dict['expression_subtype_magnoid'].append(sampleID)
        elif expr: 
            null_dict['expression_subtype_magnoid'].append(sampleID)

        expr = line[clinical_header.index('Expression_Subtype')]
        if expr == 'Squamoid':
            sets_dict['expression_subtype_squamoid'].append(sampleID)
        elif expr: 
            null_dict['expression_subtype_squamoid'].append(sampleID)        

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
            #add pts to dead at 1yrs
            if float(dtd) < 365 and dead == '1':
                sets_dict['dead_1yr'].append(sampleID)
            elif dtd:
                null_dict['dead_1yr'].append(sampleID)
                

        
        
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
        node = line[clinical_header.index('lymphnode_pathologic_spread')]
        if node == ('N1' or 'N2' or 'N3'):
            sets_dict['node_positive'].append(sampleID)
        elif node: 
            null_dict['node_positive'].append(sampleID)
        if node == ('N0'):
            sets_dict['node_negative'].append(sampleID)
        elif node: 
            null_dict['node_negative'].append(sampleID)

        
        stage = line[clinical_header.index('pathologic_stage')]    
        #add to stage >=2 set
        if stage == 'Stage IIA' or \
            stage == 'Stage IIB' or \
            stage == 'Stage IIIA' or \
            stage == 'Stage IIIB' or \
            stage == 'Stage IV':
            sets_dict['stage_gte_II'].append(sampleID)
        elif stage != '': 
            null_dict['stage_gte_II'].append(sampleID)
            
        
        #add to stage >=3 set
        if  stage == 'Stage IIIA' or \
            stage == 'Stage IIIB' or \
            stage == 'Stage IV':        
            sets_dict['stage_gte_III'].append(sampleID)
        elif stage != '': 
            null_dict['stage_gte_III'].append(sampleID)    
        
        if stage == 'Stage IV':
            sets_dict['metastatic'].append(sampleID)
        elif stage: 
            null_dict['metastatic'].append(sampleID)
            
        #check node status
        disease_status = line[clinical_header.index('primary_therapy_outcome_success')]
        if disease_status == ('Progressive Disease'):
            sets_dict['progressive_disease'].append(sampleID)
        elif disease_status: 
            null_dict['progressive_disease'].append(sampleID)
        if disease_status == ('Complete Remission/Response' or 'Partial Remission/Response'):
            sets_dict['disease_in_remission'].append(sampleID)
        elif disease_status: 
            null_dict['disease_in_remission'].append(sampleID)
        
        radiation = line[clinical_header.index('radiation_therapy')]
        if radiation == ('YES'):
            sets_dict['radiation_therapy'].append(sampleID)
        elif radiation: 
            null_dict['radiation_therapy'].append(sampleID)
        
        smoking = line[clinical_header.index('tobacco_smoking_history')]
        if smoking == ('Current reformed smoker for < or = 15 years'):
            sets_dict['reformed_smoker_lte_15yrs'].append(sampleID)
        elif smoking: 
            null_dict['reformed_smoker_lte_15yrs'].append(sampleID)
        if smoking == ('Current reformed smoker for > 15 years'):
            sets_dict['reformed_smoker_gt_15yrs'].append(sampleID)
        elif smoking: 
            null_dict['reformed_smoker_gt_15yrs'].append(sampleID)
        if smoking == ('Current smoker'):
            sets_dict['current_smoker'].append(sampleID)
        elif smoking: 
            null_dict['current_smoker'].append(sampleID)
        if smoking == ('Lifelong Non-smoker'):
            sets_dict['never_smoker'].append(sampleID)
        elif smoking: 
            null_dict['never_smoker'].append(sampleID)
        

        
        
    ca_type = 'LUAD' 
    
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