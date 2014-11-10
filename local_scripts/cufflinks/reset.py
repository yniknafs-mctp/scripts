'''
Created on Aug 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import glob
import shutil



    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("cuff_dir")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Cleaning cufflinks dir')
    cuff_dirs = glob.glob(os.path.join(args.cuff_dir, '*'))
    tot = len(cuff_dirs)
    i = 0
    for dir in cuff_dirs:
        i+=1
        if i%100 == 0:
            logging.info("Finished %d/%d" % (i, tot))
        files = glob.glob(os.path.join(dir, '*'))
        for file in files:
            base = os.path.basename(file)
            if base == 'tmp':
                shutil.rmtree(file)
            elif base == 'cufflinks.sh':
                continue
            else: 
                os.remove(file) 
                
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
