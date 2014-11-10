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
    desc_dict['age_diagnosis_lt40'].append('Pts who were younger than 40 at time of diagnosis')
    desc_dict['age_diagnosis_lt50'].append('Pts who were younger than 50 at time of diagnosis')
    desc_dict['age_diagnosis_gt65'].append('Pts who were older than 65 at time of diagnosis')
    desc_dict['age_diagnosis_gt75'].append('Pts who were older than 75 at time of diagnosis')
    desc_dict['dead_3mos'].append('Pts who died within three months of diagnosis')    
    desc_dict['dead_6mos'].append('Pts who died within six months of diagnosis')
    desc_dict['dead_9mos'].append('Pts who died within nine months of diagnosis')
    desc_dict['dead_12mos'].append('Pts who died within twelve months of diagnosis')
    desc_dict['CDE_chemo_adjuvant_alk_positive'].append('Pts who received adjuvant alk chemo')
    desc_dict['CDE_chemo_adjuvant_alk_negative'].append('Pts who did not receive adjuvant alk chemo')
    desc_dict['CDE_chemo_adjuvant_tmz_positive'].append('Pts who received adjuvant tmz chemo')
    desc_dict['CDE_chemo_adjuvant_tmz_negative'].append('Pts who did not receive adjuvant tmz chemo')
    desc_dict['CDE_chemo_alk_positive'].append('Pts who received alk chemo')
    desc_dict['CDE_chemo_alk_negative'].append('Pts who did not receive alk chemo')
    desc_dict['CDE_chemo_tmz_positive'].append('Pts who received tmz chemo')
    desc_dict['CDE_chemo_tmz_negative'].append('Pts who did not receive tmz chemo')
    desc_dict['CDE_radiation_adjuvant_positive'].append('Pts who received adjuvant radiation')
    desc_dict['CDE_radiation_adjuvant_negative'].append('Pts who did not receive adjuvant radiation')
    desc_dict['CDE_radiation_any_positive'].append('Pts who received radiation')
    desc_dict['CDE_radiation_any_negative'].append('Pts who did not receive radiation')
    desc_dict['GeneExp_classical'].append('Samples that have classical GeneExp')
    desc_dict['GeneExp_mesenchymal'].append('Samples that have mesenchymal GeneExp')
    desc_dict['GeneExp_proneural'].append('Samples that have proneural GeneExp')
    desc_dict['GeneExp_neural'].append('Samples that have neural GeneExp')
    desc_dict['Male'].append('Male patients')
    desc_dict['Female'].append('Female patients')
    desc_dict['karnofsky_lte80'].append('Pts with a Karnofsky performance score less than or equal to 80%')
    desc_dict['karnofsky_lte70'].append('Pts with a Karnofsky performance score less than or equal to 70%')
    desc_dict['node_positive'].append('Samples that are node positive')
    desc_dict['node_negative'].append('Samples that are node negative')
    desc_dict['recurrence'].append('Samples that recurrent tumors')
    desc_dict['G_CIMP_negative'].append('patients that are G-CIMP negative')
    desc_dict['G_CIMP_positive'].append('patients that are G-CIMP negative')
    
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)
        #checking age at diagnosis field 
        age_diagnosis = line[clinical_header.index('age_at_initial_pathologic_diagnosis')]
        if age_diagnosis: 
            #add to age<40 set
            if float(age_diagnosis) < 40:
                    sets_dict['age_diagnosis_lt40'].append(sampleID)
            elif age_diagnosis: 
                    null_dict['age_diagnosis_lt40'].append(sampleID)
            #add to age<50 set
            if float(age_diagnosis) < 50:
                    sets_dict['age_diagnosis_lt50'].append(sampleID)
            elif age_diagnosis: 
                    null_dict['age_diagnosis_lt50'].append(sampleID)
            #add to age<60 set
            if float(age_diagnosis) > 65:
                    sets_dict['age_diagnosis_gt65'].append(sampleID)
            elif age_diagnosis: 
                    null_dict['age_diagnosis_gt65'].append(sampleID)
            #add to age>75 set
            if float(age_diagnosis) > 75:
                    sets_dict['age_diagnosis_gt75'].append(sampleID)
            elif age_diagnosis: 
                    null_dict['age_diagnosis_gt75'].append(sampleID)
             
                
            
            
        #check CDE_chemo_adjuvant_alk field
        CDE_chemo_adjuvant_alk = line[clinical_header.index('CDE_chemo_adjuvant_alk')]
        #add pts to CDE_chemo_adjuvant_alk_positive
        if CDE_chemo_adjuvant_alk == 'TRUE': 
            sets_dict['CDE_chemo_adjuvant_alk_positive'].append(sampleID)
        elif CDE_chemo_adjuvant_alk: 
            null_dict['CDE_chemo_adjuvant_alk_positive'].append(sampleID)
        
        #add pts to CDE_chemo_adjuvant_alk_negative
        if CDE_chemo_adjuvant_alk == 'FALSE': 
            sets_dict['CDE_chemo_adjuvant_alk_negative'].append(sampleID)
        elif CDE_chemo_adjuvant_alk:
            null_dict['CDE_chemo_adjuvant_alk_negative'].append(sampleID)
        
        #check CDE_chemo_adjuvant_tmz field
        CDE_chemo_adjuvant_tmz = line[clinical_header.index('CDE_chemo_adjuvant_tmz')]
        #add pts to CDE_chemo_adjuvant_tmz_positive
        if CDE_chemo_adjuvant_tmz == 'TRUE': 
            sets_dict['CDE_chemo_adjuvant_tmz_positive'].append(sampleID)
        elif CDE_chemo_adjuvant_tmz: 
            null_dict['CDE_chemo_adjuvant_tmz_positive'].append(sampleID)
        #add pts to CDE_chemo_adjuvant_tmz_negative
        if CDE_chemo_adjuvant_tmz == 'FALSE': 
            sets_dict['CDE_chemo_adjuvant_tmz_negative'].append(sampleID)
        elif CDE_chemo_adjuvant_tmz: 
            null_dict['CDE_chemo_adjuvant_tmz_negative'].append(sampleID)
        
        #check CDE_chemo_alk field
        CDE_chemo_alk = line[clinical_header.index('CDE_chemo_alk')]
        #add pts to CDE_chemo_alk_positive
        if CDE_chemo_alk == 'TRUE': 
            sets_dict['CDE_chemo_alk_positive'].append(sampleID)
        elif CDE_chemo_alk: 
            null_dict['CDE_chemo_alk_positive'].append(sampleID)
        #add pts to CDE_chemo_alk_negative
        if CDE_chemo_alk == 'FALSE': 
            sets_dict['CDE_chemo_alk_negative'].append(sampleID)
        elif CDE_chemo_alk: 
            null_dict['CDE_chemo_alk_negative'].append(sampleID)
        
        #check CDE_chemo_tmz field
        CDE_chemo_tmz = line[clinical_header.index('CDE_chemo_tmz')]
        #add pts to CDE_chemo_tmz_positive
        if CDE_chemo_tmz == 'TRUE': 
            sets_dict['CDE_chemo_tmz_positive'].append(sampleID)
        elif CDE_chemo_tmz: 
            null_dict['CDE_chemo_tmz_positive'].append(sampleID)
        #add pts to CDE_chemo_tmz_negative
        if CDE_chemo_tmz == 'FALSE': 
            sets_dict['CDE_chemo_tmz_negative'].append(sampleID)
        elif CDE_chemo_tmz: 
            null_dict['CDE_chemo_tmz_negative'].append(sampleID)

        #check adjuvant radiation field
        CDE_radiation_adjuvant = line[clinical_header.index('CDE_radiation_adjuvant')]
        #add pts to CDE_radiation_adjuvant_positive
        if CDE_radiation_adjuvant == 'TRUE': 
            sets_dict['CDE_radiation_adjuvant_positive'].append(sampleID)
        elif CDE_radiation_adjuvant: 
            null_dict['CDE_radiation_adjuvant_positive'].append(sampleID)
        #add pts to CDE_radiation_adjuvant_negative
        if CDE_radiation_adjuvant == 'FALSE': 
            sets_dict['CDE_radiation_adjuvant_negative'].append(sampleID)
        elif CDE_radiation_adjuvant: 
            null_dict['CDE_radiation_adjuvant_negative'].append(sampleID)

        #check radiation any field
        CDE_radiation_any = line[clinical_header.index('CDE_radiation_any')]
        #add pts to CDE_radiation_any_positive
        if CDE_radiation_any == 'TRUE': 
            sets_dict['CDE_radiation_any_positive'].append(sampleID)
        elif CDE_radiation_any: 
            null_dict['CDE_radiation_any_positive'].append(sampleID)
        #add pts to CDE_chemo_tmz_negative
        if CDE_radiation_any == 'FALSE': 
            sets_dict['CDE_radiation_any_negative'].append(sampleID)
        elif CDE_radiation_any: 
            null_dict['CDE_radiation_any_negative'].append(sampleID)

        
        
        
        #check Days to date of death field
        dtd = line[clinical_header.index('_OS')]
        dead = line[clinical_header.index('_EVENT')]
        if dtd: 
            #add pts to dead at 3mos
            if float(dtd) < 90 and dead == '1': 
                sets_dict['dead_3mos'].append(sampleID)
            elif dtd: 
                null_dict['dead_3mos'].append(sampleID)
            #add pts to dead at 6mos
            if float(dtd) < 180 and dead == '1':
                sets_dict['dead_6mos'].append(sampleID)
            elif dtd: 
                null_dict['dead_6mos'].append(sampleID)
            #add pts to dead at 9mos
            if float(dtd) < 270 and dead == '1':
                sets_dict['dead_9mos'].append(sampleID)
            elif dtd:
                null_dict['dead_9mos'].append(sampleID)
            #add pts to dead at 12mos
            if float(dtd) < 360 and dead == '1':
                sets_dict['dead_12mos'].append(sampleID)
            elif dtd:
                null_dict['dead_12mos'].append(sampleID)
      
        
        #check G_CIMP status
        G_CIMP = line[clinical_header.index('G_CIMP_STATUS')]
        #add pts to G_CIMP_positive
        if G_CIMP == 'G-CIMP': 
            sets_dict['G_CIMP_positive'].append(sampleID)
        elif G_CIMP: 
            null_dict['G_CIMP_positive'].append(sampleID)
        #add pts to G_CIMP_negative
        if G_CIMP == 'NON G-CIMP': 
            sets_dict['G_CIMP_negative'].append(sampleID)
        elif G_CIMP: 
            null_dict['G_CIMP_negative'].append(sampleID)

        #check GeneExp_Subtype
        GeneExp = line[clinical_header.index('GeneExp_Subtype')]
        #add pts to GeneExp_classical
        if GeneExp == 'Classical': 
            sets_dict['GeneExp_classical'].append(sampleID)
        elif GeneExp: 
            null_dict['GeneExp_classical'].append(sampleID)
        #add pts to GeneExp_mesenchymal
        if GeneExp == 'Mesenchymal': 
            sets_dict['GeneExp_mesenchymal'].append(sampleID)
        elif GeneExp: 
            null_dict['GeneExp_mesenchymal'].append(sampleID)
        #add pts to GeneExp_proneural
        if GeneExp == 'Proneural': 
            sets_dict['GeneExp_proneural'].append(sampleID)
        elif GeneExp: 
            null_dict['GeneExp_proneural'].append(sampleID)
        #add pts to GeneExp_neural
        if GeneExp == 'Neural': 
            sets_dict['GeneExp_neural'].append(sampleID)
        elif GeneExp: 
            null_dict['GeneExp_neural'].append(sampleID)

        #check gender
        gender = line[clinical_header.index('gender')]
        #add pts to male
        if gender == 'MALE': 
            sets_dict['Male'].append(sampleID)
        elif gender: 
            null_dict['Male'].append(sampleID)
        #add pts to female
        if gender == 'FEMALE': 
            sets_dict['Female'].append(sampleID)
        elif gender: 
            null_dict['Female'].append(sampleID)
            
        #check karnofsky score
        karnofsky = line[clinical_header.index('karnofsky_performance_score')]
        if karnofsky:
            #add pts to karnofsky_lte80
            if float(karnofsky) <= 80: 
                sets_dict['karnofsky_lte80'].append(sampleID)
            elif karnofsky:
                null_dict['karnofsky_lte80'].append(sampleID)
            #add pts to karnofsky_lte70
            if float(karnofsky) <= 70: 
                sets_dict['karnofsky_lte70'].append(sampleID)
            elif karnofsky:
                null_dict['karnofsky_lte70'].append(sampleID)
        
        #check recurrence
        recurrence = line[clinical_header.index('sample_type')]
        #add pts to male
        if recurrence == 'Recurrent Tumor': 
            sets_dict['recurrence'].append(sampleID)
        elif recurrence: 
            null_dict['recurrence'].append(sampleID)
    
    smt_file = open('GBM_clinical_sets.smt', 'w')
    header_row = ['Set name', 'Set description'] + universe
    smt_file.write('\t'.join(map(str, header_row)))
    smt_file.write('\n')
    for key in sets_dict.iterkeys():
        desc = 'GBM ' +desc_dict[key][0]
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
        name = 'GBM ' + key.replace('_', ' ')
        line = [name] + [desc] + membership
        smt_file.write('\t'.join(map(str,line)))
        smt_file.write('\n')
    smt_file.close()
        

if __name__ == '__main__':
    sys.exit(main())