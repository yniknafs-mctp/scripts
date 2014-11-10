'''
Created on Aug 4, 2014
@author yniknafs
'''


import os
import sys
import logging
import subprocess
import re


'''
run on hpc cluster to check the number of free nodes 
'''
    
def main():
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    tot_nodes = 384
    used_nodes = 0
    qstat = subprocess.check_output('qstat -R| tail -n +6', shell=True).strip().split('\n')
    for line in qstat: 
        fields = re.split('\s+', line)
        if fields[7] == 'R':
            used_nodes+=int(fields[4])
    
    
    logging.info('LOAD:\t %s/%s' % (used_nodes,tot_nodes))
    logging.info('FREE:\t %d' % (tot_nodes - used_nodes))
    
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())

