'''
Created on October 28, 2013
@author: yniknafs    
'''
import os
import sys
import argparse
import logging
import numpy as np
import csv


def CNV_pts(headers, line, last_pheno_col, amp_thresh, del_thresh):
    
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    #logging.info('Identifying CNV samples')
    
    
    cnv_pts = []
    null_pts = []
    
    peak_type = line[0].split('Peak')[0]

    if peak_type == "Amplification":
        thresh = float(amp_thresh)
        multiplier = 1
    elif peak_type == "Deletion":
        thresh = float(del_thresh)
        multiplier = -1
    else:
        logging.debug("Not a valid peak type in CNV_sample identification")
        sys.exit()
    
    #categorize the samples based on whether or not they are amplified
    for x in xrange(last_pheno_col, len(line)):
        #use this product so that you don't have to have multiple if statements for directionality
        product = float(line[x])*multiplier
        #print product
        if product > thresh:
            if headers[x] != []:
                cnv_pts.append(headers[x])
        else:
            if headers[x] != []:
                null_pts.append(headers[x])
    return cnv_pts, null_pts
                
    
    
    
    
    
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
