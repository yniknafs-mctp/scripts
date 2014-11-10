'''
Created on Aug 21, 2014
@author yniknafs
'''


import os
import sys
import argparse
import logging
import paramiko
import subprocess
import re



    
def main():
    # parse command line
    parser = argparse.ArgumentParser()
    
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    
    tot_nodes_all = 0
    #check hpc
    tot_nodes = 384
    used_nodes = 0
    HOST="mctp-hpc-login3"
    COMMAND="qstat -R| tail -n +6"
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    
    qstat = ssh.stdout.readlines()
    for line in qstat: 
        fields = re.split('\s+', line)
        if fields[7] == 'R':
            used_nodes+=int(fields[4])
    logging.info('-------HPC------')
    logging.info('LOAD:\t %s/%s' % (used_nodes,tot_nodes))
    logging.info('FREE:\t %d' % (tot_nodes - used_nodes))
    tot_nodes_all += tot_nodes-used_nodes
    #check pb8
    tot_nodes = 64
    used_nodes = 0
    HOST="pathbio-8"
    COMMAND="ps aux| tail -n +2"
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    for line in ps: 
        fields = re.split('\s+', line)
        if float(fields[2]) > 40:
            used_nodes+=1
    logging.info('-------PB8------')
    logging.info('LOAD:\t %s/%s' % (used_nodes,tot_nodes))
    logging.info('FREE:\t %d' % (tot_nodes - used_nodes))
    tot_nodes_all += tot_nodes-used_nodes
    #check pb9
    tot_nodes = 64
    used_nodes = 0
    HOST="pathbio-9"
    COMMAND="ps aux| tail -n +2"
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    for line in ps: 
        fields = re.split('\s+', line)
        if float(fields[2]) > 40:
            used_nodes+=1
    logging.info('-------PB9------')
    logging.info('LOAD:\t %s/%s' % (used_nodes,tot_nodes))
    logging.info('FREE:\t %d' % (tot_nodes - used_nodes))
    tot_nodes_all += tot_nodes-used_nodes
    logging.info('----------------')
    logging.info('Total free nodes: %d' % tot_nodes_all)
    
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())
