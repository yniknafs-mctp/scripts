'''
Created on Aug 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import re
import subprocess




    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("user")
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    
    logging.info('Killing all jobs for %s' %args.user)
    qstat = subprocess.check_output('qstat -u yniknafs -R| tail -n +6', shell=True).strip().split('\n')
    for line in qstat: 
        fields = re.split('\s+', line)
        jobid = fields[0].split('.')[0]
        arg = 'qdel %s' % jobid
        logging.debug(arg)
        subprocess.call(arg, shell=True)
       
    return 0
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
