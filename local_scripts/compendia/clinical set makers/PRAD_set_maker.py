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
    desc_dict['gleason_7'].append('Patients with overall gleason score of 7')
    desc_dict['gleason_gte_8'].append('Patients with overall gleason score greater than 8')
    desc_dict['gleason_gte_9'].append('Patients with overall gleason score greater than 9')
    desc_dict['psa_lt_10'].append('Patients with PSA less than 10 at time of sample collection')
    desc_dict['psa_gte_10'].append('Patients with PSA greater than 10 at time of sample collection')
    desc_dict['psa_gte_20'].append('Patients with PSA greater than 20 at time of sample collection')
    desc_dict['node_positive'].append('Pts with node positive disease')
    desc_dict['node_negative'].append('Pts with node negative disease')

            
        
    #loop through clincal file, adding samples to various sets
    for line in clinical_file:
        line = line.strip().split('\t')
        #scan each line for elements used to make sets and add samples to sets
        
        sampleID = line[clinical_header.index('sampleID')]
        universe.append(sampleID)

        
        gleason = line[clinical_header.index('gleason_score')]
        #add to stage >=2 set
        if gleason == '7':
            sets_dict['gleason_7'].append(sampleID)
        elif gleason != '': 
            null_dict['gleason_7'].append(sampleID)
        if gleason == '9' or gleason == '10' or gleason == '8':
            sets_dict['gleason_gte_8'].append(sampleID)
        elif gleason != '': 
            null_dict['gleason_gte_8'].append(sampleID)
        if gleason == '9' or gleason == '10':
            sets_dict['gleason_gte_9'].append(sampleID)
        elif gleason != '': 
            null_dict['gleason_gte_9'].append(sampleID)
        
        psa = line[clinical_header.index('psa_result_most_recent')]
        if psa:
            psa = float(psa)
            if psa < 10:
                sets_dict['psa_lt_10'].append(sampleID)
            elif psa: 
                null_dict['psa_lt_10'].append(sampleID)
            if psa >= 10:
                sets_dict['psa_gte_10'].append(sampleID)
            elif psa: 
                null_dict['psa_gte_10'].append(sampleID)
            if psa > 20:
                sets_dict['psa_gte_20'].append(sampleID)
            elif psa: 
                null_dict['psa_gte_20'].append(sampleID)    
        
        
        #check node status
        node = line[clinical_header.index('pathologic_N')]
        if node == 'N1':
            sets_dict['node_positive'].append(sampleID)
        elif node: 
            null_dict['node_positive'].append(sampleID)
        if node == 'N0':
            sets_dict['node_negative'].append(sampleID)
        elif node: 
            null_dict['node_negative'].append(sampleID)
        

        

  
  
    ca_type = 'PRAD' 
    
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