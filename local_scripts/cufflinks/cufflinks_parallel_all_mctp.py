    
'''
Created on May 10, 2014
@author yniknafs
'''

import os
import sys
import argparse
import logging
import subprocess
import glob
import time
import re

COMPUTE_GENE_EXPR_SCRIPT = '/mctp/users/yniknafs/scripts/workspace/assemblyline/assemblyline/pipeline/compute_gene_expression.py'
USER = 'yniknafs'

def hpc_load():
    #check hpc
    HOST="mctp-hpc-login3"
    COMMAND='qstat -u %s| grep " R " | wc -l' % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    used_nodes = 0
    qstat = ssh.stdout.read()
#     for line in qstat: 
#         fields = re.split('\s+', line)
#         if fields[7] == 'R':
#             used_nodes+=int(fields[4])
    return qstat

def hpc_job(script):
    #kick off job on hpc
    HOST="mctp-hpc-login3"
    COMMAND='qsub %s' % script
    subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)

def pb8_load():
    #check pb8
    used_nodes =0
    HOST="pathbio-8"
    COMMAND="ps aux|grep %s" % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    for line in ps: 
        fields = re.split('\s+', line)
        if float(fields[2]) > 40:
            used_nodes+=1
    return used_nodes

def pb8_job(script):
    #kick off job on hpc
    HOST="pathbio-8"
    COMMAND='sh %s' % script
    subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)


def pb9_load():
    #check pb9
    used_nodes =0
    HOST="pathbio-8"
    COMMAND="ps aux|grep %s" % USER
    ssh = subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    ps = ssh.stdout.readlines()
    for line in ps: 
        fields = re.split('\s+', line)
        if float(fields[2]) > 40:
            used_nodes+=1
    return used_nodes

def pb9_job(script):
    #kick off job on hpc
    HOST="pathbio-9"
    COMMAND='sh %s' % script
    subprocess.Popen(["ssh", "%s" % HOST, COMMAND],
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)



def sleeper(hpc_procs, pb8_procs, pb9_procs):
    while hpc_load() >= hpc_procs and pb8_procs >= pb8_load() and pb9_procs > pb9_load(): 
        time.sleep(30)
    free_agents = []
    if hpc_load() < hpc_procs: 
        free_agents.append('hpc')
    if pb8_load() < pb8_procs: 
        free_agents.append('pb8')
    if pb9_load() < pb9_procs: 
        free_agents.append('pb9')
    logging.debug(free_agents)
    logging.debug(hpc_load())
    logging.debug(hpc_procs)
    return free_agents
    

def main():
    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("config_xml")
    parser.add_argument("library_table")
    parser.add_argument("gtf_file")
    parser.add_argument("out_dir")
    parser.add_argument("-hpc", dest = 'hpc',
                    default = 4,
                    help = 'number of HPC processors to use')
    parser.add_argument("-pb8", dest = 'pb8',
                    default = 4,
                    help = 'number of pb8 processors to use')
    parser.add_argument("-pb9", dest = 'pb9',
                    default = 4,
                    help = 'number of pb9 processors to use')
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG,
                      format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
     
    args.hpc = int(args.hpc)
    args.pb8 = int(args.pb8)
    args.pb9 = int(args.pb9)
    #run the script to make the PBS scripts
    script_gen_args = [
                       'python',
                       COMPUTE_GENE_EXPR_SCRIPT,
                       '-o',
                       args.out_dir, 
                       args.config_xml,
                       args.library_table,
                       args.gtf_file,
                       '--mode',
                       'cufflinks' 
                       ]
      
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
      
#     logging.debug('Generating PBS scripts...')
#     subprocess.call(script_gen_args)
      
    sample_dirs = glob.glob(os.path.join(args.out_dir, '*'))
      
    #kick off jobs on cluster run cufflinks
    logging.debug('Running %d jobs total...' % len(sample_dirs))
    logging.debug('Using %d cores on HPC' % args.hpc)
    logging.debug('Using %d cores on PB8' % args.pb8)
    logging.debug('Using %d cores on PB9' % args.pb9)
    sleeper(args.hpc, args.pb8, args.pb9)
#     i=0
#     for sample in sample_dirs:
#         i+=1
#         script = os.path.join(sample, 'cufflinks.sh')
#         free_agents = sleeper(args.hpc, args.pb8, args.pb9)
#         if 'hpc' in free_agents:
#             logging.debug('Running job %d/%d on HPC: %s' % (i, len(sample_dirs), os.path.basename(sample)))
#             hpc_job(script)
#         elif 'pb8' in free_agents:
#             logging.debug('Running job %d/%d on PB8: %s' % (i, len(sample_dirs), os.path.basename(sample)))
#             pb8_job(script)
#         elif 'pb9' in free_agents:
#             logging.debug('Running job %d/%d on PB9: %s' % (i, len(sample_dirs), os.path.basename(sample)))
#             pb9_job(script)
            
    
    return 0

if __name__ == '__main__': 
    sys.exit(main())    

